
import os
import shutil

class handle_workdir:

    def __init__( self, dp ):

        # dp := script's directory path 
        
        if dp[-1] != '/': dp+='/'
        self.dp         = dp

        self.wd_prepend = 'siracha_workdir_'

        # Assign working directory name
        self.find_workdirs()
        self.wd_name+='/'

        # Make working directory
        os.mkdir( self.dp+self.wd_name )
        self.wd_fullpath = self.dp+self.wd_name

        # Move rrtmg files
        shutil.copy( 
                  "/home/edward/lib/siracha/rrtmg_sw"
                , self.dp+self.wd_name+"rrtmg_sw"
                )
        shutil.copy( 
                  "/home/edward/lib/siracha/rrtmg_lw"
                , self.dp+self.wd_name+"rrtmg_lw"
                )

        # Now the workdir is ready. Clean up with handle_workdir.clean()

    def get_full_path( self ):

        return self.wd_fullpath

    def clean( self ):

        shutil.rmtree( self.dp+self.wd_name )

    def find_workdirs( self ):

        workdir_names = []
        workdir_nums  = []

        dp_contents = os.listdir( self.dp )

        for f in dp_contents:
            if f.find( self.wd_prepend )>=0: workdir_names.append(f)

        if len(workdir_names)==0:
            self.wd_name = self.wd_prepend+'1'
            return

        # Now we want the numbers of the existing working directories

        for w in workdir_names:
            workdir_nums.append( int(w[len(self.wd_prepend):]) )

        workdir_nums.sort()

        if len(workdir_nums)==max(workdir_nums):
            self.wd_name    = self.wd_prepend+str(max(workdir_nums)+1)
            return
        else:
            for n in range(1,max(workdir_nums)):
                if n not in workdir_nums: 
                    self.wd_name = self.wd_prepend+str(n)
                    return






