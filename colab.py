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

        subprocess.run('git clone https://github.com/michaelJwilson/desihigh.git', shell=True, check=True)    
        # subprocess.run('pip install -r desihigh/requirements.txt', shell=True, check=True)

        try:
            import desihigh

        except:
            emessage = 'Failed to setup colab.  Create a ticket at https://github.com/michaelJwilson/desihigh.git.'

            raise RuntimeError(emessage)


def save_colab():
    if 'COLAB_GPU' in os.environ:
        # drive.mount(drive, force_remount=True)

        drive.flush_and_unmount()
