import numpy as np
import modflowapi


def get_model_variables(modelname):
    model_variables = [
    # gwf
    ("insto", modelname),
    ("iss", modelname),
    ("x", modelname),
    ("xold", modelname),
    ("rhs", modelname),

    # dis
    ("area", modelname, "dis"),
    ("top", modelname, "dis"),
    ("bot", modelname, "dis"),

    # npf
    ("sat", modelname, "npf"),
    ("condsat", modelname, "npf"),
    ("derv_mult", modelname, "npf"),

    # swi
    ("zeta", modelname, "swi"),
    ("hcof", modelname, "swi"),
    ("rhs", modelname, "swi"),
    ]
    return model_variables


class SwiAPI:

    def __init__(self, libmf6, ws, modelname, salt_modelname=None):
        self.mf6api =  modflowapi.ModflowApi(libmf6, working_directory=ws)
        self.modelname = modelname
        self.salt_modelname = salt_modelname
        self.initialize()
        self.sim_pointer = {}
        self.model_fresh_pointer = {}
        self.model_salt_pointer = {}
        self.zeta_all = []
        return
    
    def initialize(self):
        self.rhof = 1000.
        self.rhos = 1025.
        drho = self.rhos - self.rhof
        self.alpha_f = self.rhof / drho
        self.alpha_s = self.rhos / drho
        self.head_saltwater = 0.
    
    def create_pointers(self, verbose):
        """Create pointers to modflow variables and store in
           self.sim_pointer, self.model_fresh_pointer, 
           self.model_salt_pointer.
        """
        sim_variables = [

            # tdis
            ("delt", 'tdis'),
            ("nper", 'tdis'),

            # sln
            ("ia", "sln_1"),
            ("ja", "sln_1"),
            ("amat", "sln_1"),
        ]

        # add list of simulation pointers to self.sim_pointer dict
        for v in sim_variables:
            self.add_pointer(self.sim_pointer, v, verbose)

        # add list of pointers to fresh model dictionary
        for v in get_model_variables(self.modelname):
            self.add_pointer(self.model_fresh_pointer, v, verbose)

        # add list of pointers to salt model dictionary
        if self.salt_modelname is not None:
            for v in get_model_variables(self.salt_modelname):
                self.add_pointer(self.model_salt_pointer, v, verbose)

        # check for storage and add storage pointers
        if self.model_fresh_pointer["insto"] > 0:
            v = ("sy", self.modelname, "sto")
            self.add_pointer(self.model_fresh_pointer, v, verbose)

        # add list of pointers to salt model dictionary
        if self.salt_modelname is not None:
            if self.model_salt_pointer["insto"] > 0:
                v = ("sy", self.salt_modelname, "sto")
                self.add_pointer(self.model_salt_pointer, v, verbose)

        # set derivative multiplier for freshwater
        derv_mult = self.model_fresh_pointer["derv_mult"]
        derv_mult[:] = self.alpha_f + 1.0

        # set derivative multiplier for saltwater
        if self.salt_modelname is not None:
            derv_mult = self.model_salt_pointer["derv_mult"]
            derv_mult[:] = self.alpha_s

    def add_pointer(self, pointer_dict, variable, verbose):
        variable_name = variable[0]
        args = [item.upper() for item in variable]
        tag = self.mf6api.get_var_address(*args)
        if verbose:
            print(f"Accessing pointer using tag: {tag}")
        pointer = self.mf6api.get_value_ptr(tag)
        pointer_dict[variable_name] = pointer

    def print_pointers(self):
        """Print pointers to MODFLOW variables"""
        for pointer_dict in [
            self.sim_pointer, 
            self.model_fresh_pointer, 
            self.model_salt_pointer]:
            for variable_name in pointer_dict:
                pointer = pointer_dict[variable_name]
                print(f"{variable_name}: {pointer}")

    def update_zeta(self):
        """Recalculate zeta using the pointer to the modflow 
        freshwater head and saltwater head (if saltwater model
        is present)."""
        zeta_fresh = self.model_fresh_pointer["zeta"]
        self.zeta_last = zeta_fresh.copy()
        head_freshwater = self.model_fresh_pointer["x"]
        head_saltwater = self.head_saltwater
        if self.salt_modelname is not None:
            head_saltwater = self.model_salt_pointer["x"]
            zeta_salt = self.model_salt_pointer["zeta"]
        zeta = self.alpha_s * head_saltwater - self.alpha_f * head_freshwater
        zeta_fresh[:] = zeta
        if self.salt_modelname is not None:
            zeta_salt[:] = zeta
        print (f"Setting {zeta=}")

    def formulate(self, dt):
        is_transient = self.model_fresh_pointer["iss"] == 0
        if is_transient:
            self.formulate_storage(dt, self.model_fresh_pointer, self.alpha_f)
            if self.salt_modelname is not None:
                # flip sign for saltwater model
                self.formulate_storage(dt, self.model_salt_pointer, -self.alpha_s)
        return

    def formulate_storage(self, dt, model_pointer, alpha):
        """
        Calculate hcof and rhs in the modflow swi package to 
        account for the change in freshwater storage.

            hcof = - alpha_f * area * S_zeta / delt
            rhs = - alpha_f * area * S_zeta / delt * h_old 

        This function is intended to work for either the freshwater
        model or the saltwater model (if present)
        """
        head = model_pointer["x"]
        hold = model_pointer["xold"]
        hcof = model_pointer["hcof"]
        rhs = model_pointer["rhs"]
        zeta = model_pointer["zeta"]
        top = model_pointer["top"]
        bot = model_pointer["bot"]
        area = model_pointer["area"]
        sy = model_pointer["sy"]
        s_zeta = np.zeros(head.shape)
        idx = (zeta > bot) & (zeta <= top)
        s_zeta[idx] = sy[idx]
        sc = - alpha / dt * np.multiply(area, s_zeta)
        hcof[:] = sc
        rhs[:] = sc * hold

    def run(self, maxiter, verbose=True):
        if verbose:
            print("Initializing mf6...")

        self.mf6api.initialize()
        self.create_pointers(verbose)
        # self.print_pointers()

        current_time = 0.
        end_time = self.mf6api.get_end_time()
        if verbose:
            print(f"Simulation end time = {end_time}...")

        # update zeta using the initial head
        self.update_zeta()

        while current_time < end_time:
            
            if verbose:
                print(f"\n  Solving for time {current_time}")
            # dt = self.mf6api.get_time_step()
            dt = self.sim_pointer["delt"]

            if verbose:
                print(f"  Prepare time step with dt={dt}...")
            self.mf6api.prepare_time_step(dt)
            
            kiter = 0
            if verbose:
                print("  Prepare solve...")
            self.mf6api.prepare_solve(1)

            while kiter < maxiter:

                if verbose:
                    print(f"    Solve...(kiter={kiter})")
                self.formulate(dt)
                has_converged = self.mf6api.solve(1)

                # update zeta using the recent head iterate
                self.update_zeta()

                if has_converged:
                    break
                kiter += 1

                # swiapi.print_pointers()

            # save zeta for this time step
            self.zeta_all.append(self.zeta_last)

            if verbose:
                print("  Finalize solve...")
            self.mf6api.finalize_solve(1)
            
            if verbose:
                print("  Finalize time step...")
            self.mf6api.finalize_time_step()
            current_time = self.mf6api.get_current_time()

        if verbose:
            print ("Finalizing mf6...")
        self.mf6api.finalize()