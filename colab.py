import sys
import subprocess

from   google.colab import drive


def init_colab():
    drive='/content/drive/MyDrive/'

    drive.mount(drive, force_remount=True)

    drive = drive + '/MyDrive/'

    os.chdir(drive)

    subprocess.run('git clone https://github.com/michaelJwilson/desihigh.git', shell=True, check=True)    
    subprocess.run('pip install -r desihigh/requirements.txt', shell=True, check=True)

    sys.path.append(drive + '/desihigh/')

def save_colab():
    drive='/content/drive/MyDrive/'

    drive.mount(drive, force_remount=True)

    drive.flush_and_unmount()
