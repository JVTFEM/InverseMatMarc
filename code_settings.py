#This file allows for ease of Marc version selection, data selection, etc
class pipeline_files:
    def __init__(self):
        #main files
        mentat_version  = "2020"
        exp_model       = "sym_test_exp"
        numerical_model = "sym_test_fem"
        sub_folder      = "sph_mid"
        track_surface   = "skin"
        sample          = "sample"
        interp_direc    = 0 #0 for fem to exp, 1 for exp to fem
        interp_method   = "cubic"
        
        #Method for interpolating in the radial bassi function
        self.interp_direc  = interp_direc
        self.interp_method = interp_method
        #Name of sample evaluated
        self.sample = sample
        #name of surface for node tracking
        self.exp_track_surface = track_surface
        self.nm_track_surface  = track_surface
        #name of track surface nodes
        self.exp_track_surface_nodes = track_surface + '_nodes'
        self.nm_track_surface_nodes  = track_surface + '_nodes'
        #Assign the version of mentat to perform analysis
        self.mentat_version = mentat_version
        #Subfolder containing all the required files
        self.sub_folder = sub_folder
        #Add the neccessary tags to the experimental model files
        self.exp_mud = exp_model + '.mud'
        self.exp_t16 = exp_model + '_job1.t16'
        self.exp_dat = exp_model + '_job1.dat'
        self.exp_sts = exp_model + '_job1.sts'
        #Add the neccessary tags to the numerical model files
        self.nm_mud = numerical_model + '.mud'
        self.nm_t16 = numerical_model + '_job1.t16'
        self.nm_dat = numerical_model + '_job1.dat'
        self.nm_sts = numerical_model + '_job1.sts'

    def crnt_fldr(self,bf):   
        import pathlib
        self.curnt_folder_bck_slsh = str(pathlib.Path.cwd())
        self.curnt_folder_fwd_slsh = str(pathlib.Path.cwd()).replace("\\","/")
        if bf == "bwd":
            return self.curnt_folder_bck_slsh
        else:
            return self.curnt_folder_fwd_slsh
    