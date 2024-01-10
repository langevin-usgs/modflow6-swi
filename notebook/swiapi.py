import numpy as np
import modflowapi


class SwiAPi:

    def __init__(self, libmf6, ws, modelname):
        self.mf6api =  modflowapi.ModflowApi(libmf6, working_directory=ws)
        self.modelname = modelname
        self.initialize()
        self.api_pointer = {}
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
           self.api_pointer.
        """
        api_pointers = [

            # tdis
            ("delt", 'tdis'),
            ("nper", 'tdis'),

            # gwf
            ("insto", self.modelname),
            ("iss", self.modelname),
            ("x", self.modelname),
            ("xold", self.modelname),
            ("rhs", self.modelname),

            # sln
            ("ia", "sln_1"),
            ("ja", "sln_1"),
            ("amat", "sln_1"),

            # dis
            ("area", self.modelname, "dis"),
            ("top", self.modelname, "dis"),
            ("bot", self.modelname, "dis"),

            # npf
            ("sat", self.modelname, "npf"),
            ("condsat", self.modelname, "npf"),

            # swi
            ("zeta", self.modelname, "swi"),
            ("hcof", self.modelname, "swi"),
            ("rhs", self.modelname, "swi"),
        ]

        # add list of pointers to self.api_pointer dict
        for api_pointer in api_pointers:
            self.add_pointer(api_pointer, verbose)

        # check for storage and add storage pointers
        if self.api_pointer["insto"] > 0:
            api_pointer = ("sy", self.modelname, "sto")
            self.add_pointer(api_pointer, verbose)

    def add_pointer(self, api_pointer, verbose):
        variable_name = api_pointer[0]
        args = [item.upper() for item in api_pointer]
        tag = self.mf6api.get_var_address(*args)
        if verbose:
            print(f"Accessing pointer using tag: {tag}")
        pointer = self.mf6api.get_value_ptr(tag)
        self.api_pointer[variable_name] = pointer

    def print_pointers(self):
        """Print pointers stored in self.api_pointer"""
        for variable_name in self.api_pointer:
            pointer = self.api_pointer[variable_name]
            print(f"{variable_name}: {pointer}")

    def update_zeta(self):
        """Recalculate zeta using the pointer to the modflow head"""
        zeta = self.api_pointer["zeta"]
        self.zeta_last = zeta.copy()
        head_freshwater = self.api_pointer["x"]
        zeta[:] = (
            self.alpha_s * self.head_saltwater - self.alpha_f * head_freshwater
        )

    def formulate(self, dt):
        is_transient = self.api_pointer["iss"] == 0
        if is_transient:
            self.formulate_storage(dt)
        return

    def formulate_storage(self, dt):
        """
        Calculate hcof and rhs in the modflow swi package to 
        account for the change in freshwater storage.

            hcof = - alpha_f * area * S_zeta / delt
            rhs = - alpha_f * area * S_zeta / delt * h_old 

        """
        head = self.api_pointer["x"]
        hold = self.api_pointer["xold"]
        hcof = self.api_pointer["hcof"]
        rhs = self.api_pointer["rhs"]
        zeta = self.api_pointer["zeta"]
        top = self.api_pointer["top"]
        bot = self.api_pointer["bot"]
        area = self.api_pointer["area"]
        sy = self.api_pointer["sy"]
        s_zeta = np.zeros(head.shape)
        idx = (zeta > bot) & (zeta <= top)
        s_zeta[idx] = sy[idx]
        sc = - self.alpha_f / dt * np.multiply(area, s_zeta)
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
            dt = self.api_pointer["delt"]

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