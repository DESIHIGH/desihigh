import os
import sys
import subprocess


drivepath='/content/drive/'

if 'COLAB_GPU' in os.environ:
    from   google.colab import drive

    
    # 'Hmmm, seems you\'re not in colab :)  Try again later.'
    drive.mount(drivepath, force_remount=True)

    mydrive = drivepath + '/MyDrive/'
    
    sys.path.append(mydrive + '/desihigh/')
    os.chdir(mydrive)
    
    try:
        import  desihigh
        
    except:
        print('Failed to import desihigh; Cloning.')

        subprocess.run('git clone https://github.com/michaelJwilson/desihigh.git --depth=1', shell=True, check=True)    

        try:
            sys.path.append(mydrive + '/desihigh/')
            
            import  desihigh

            print('Successfully cloned DESI High to Google Drive.')

        except Exception as EE:
            emessage = f'Failed to setup DESI High @ colab.  Please create a ticket at https://github.com/michaelJwilson/desihigh.git and include:\n\n{EE}'

            raise  RuntimeError(emessage)

else:
    print('It appears you are not on Google Colab (!)')
