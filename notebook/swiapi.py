import numpy as np


class SwiAPi:

    def __init__(self, mfapi, modelname):
        self.mfapi = mfapi
        self.modelname = modelname
        self.initialize()
        self.api_pointer = {}
        return
    
    def initialize(self):
        self.rhof = 1000.
        self.rhos = 1025.
        drho = self.rhos - self.rhof
        self.alpha_f = self.rhof / drho
        self.alpha_s = self.rhos / drho
        self.head_saltwater = 0.
    
    def create_pointers(self):
        """Create pointers to modflow variables and store in
           self.api_pointer.
        """
        api_pointers = [

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
            self.add_pointer(api_pointer)

        # check for storage and add storage pointers
        if self.api_pointer["insto"] > 0:
            api_pointer = ("sy", self.modelname, "sto")
            self.add_pointer(api_pointer)

    def add_pointer(self, api_pointer):
        variable_name = api_pointer[0]
        args = [item.upper() for item in api_pointer]
        tag = self.mfapi.get_var_address(*args)
        print(f"Accessing pointer using tag: {tag}")
        pointer = self.mfapi.get_value_ptr(tag)
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
        sc = - self.alpha_f * area * s_zeta / dt
        hcof[:] = sc
        rhs[:] = sc * hold
