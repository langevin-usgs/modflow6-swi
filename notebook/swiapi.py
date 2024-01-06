

class SwiAPi:

    def __init__(self, mfapi, modelname):
        self.mfapi = mfapi
        self.modelname = modelname
        self.initialize()
        return
    
    def initialize(self):
        self.rhof = 1000.
        self.rhos = 1025.
        drho = self.rhos - self.rhof
        self.alpha_f = self.rhof / drho
        self.alpha_s = self.rhos / drho
        self.head_saltwater = 0.
    
    def create_pointers(self):
        """Create pointers to modflow variables"""
        api_pointers = [

            # gwf
            ("x", self.modelname),
            ("rhs", self.modelname),

            # sln
            ("ia", "sln_1"),
            ("ja", "sln_1"),
            ("amat", "sln_1"),

            # npf
            ("sat", self.modelname, "npf"),
            ("zeta", self.modelname, "npf"),
            ("condsat", self.modelname, "npf"),
        ]
        self.api_pointer = {}
        for api_pointer in api_pointers:
            variable_name = api_pointer[0]
            args = [item.upper() for item in api_pointer]
            tag = self.mfapi.get_var_address(*args)
            print(f"Accessing pointer using tag: {tag}")
            pointer = self.mfapi.get_value_ptr(tag)
            self.api_pointer[variable_name] = pointer

    def print_pointers(self):
        for variable_name in self.api_pointer:
            pointer = self.api_pointer[variable_name]
            print(f"{variable_name}: {pointer}")

    def update_zeta(self):
        zeta = self.api_pointer["zeta"]
        self.zeta_last = zeta.copy()
        head_freshwater = self.api_pointer["x"]
        zeta[:] = (
            self.alpha_s * self.head_saltwater - self.alpha_f * head_freshwater
        )

