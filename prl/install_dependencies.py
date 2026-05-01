import sys
import subprocess
import os
import platform
import urllib.request
import zipfile
import importlib

def install_python_packages():
    packages = ["numpy", "matplotlib", "scipy", "pdf2image"]
    for package in packages:
        try:
            importlib.import_module(package)
            print(f"{package} is already installed.")
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def check_poppler():
    try:
        # Check if pdftoppm is accessible
        subprocess.check_output(["pdftoppm", "-h"], stderr=subprocess.STDOUT)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def install_poppler():
    print("Poppler is missing. Attempting to install...")
    os_name = platform.system()
    
    if os_name == "Linux":
        print("Running: sudo apt-get install -y poppler-utils")
        try:
            subprocess.check_call(["sudo", "apt-get", "install", "-y", "poppler-utils"])
            print("Poppler installed successfully.")
        except subprocess.CalledProcessError:
            print("Failed to install Poppler. You may need to run the script with appropriate permissions.")
            
    elif os_name == "Darwin":
        print("Running: brew install poppler")
        try:
            subprocess.check_call(["brew", "install", "poppler"])
            print("Poppler installed successfully.")
        except subprocess.CalledProcessError:
            print("Failed to install Poppler. Please make sure Homebrew is installed.")
            
    elif os_name == "Windows":
        poppler_url = "https://github.com/oschwartz10612/poppler-windows/releases/download/v24.02.0-0/Release-24.02.0-0.zip"
        poppler_zip = "poppler.zip"
        poppler_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "poppler")
        
        if not os.path.exists(poppler_dir):
            print(f"Downloading Poppler from {poppler_url}...")
            urllib.request.urlretrieve(poppler_url, poppler_zip)
            
            print("Extracting Poppler...")
            with zipfile.ZipFile(poppler_zip, 'r') as zip_ref:
                zip_ref.extractall(poppler_dir)
            os.remove(poppler_zip)
            
        # Find the bin directory inside extracted poppler
        bin_dir = None
        for root, dirs, files in os.walk(poppler_dir):
            if "pdftoppm.exe" in files:
                bin_dir = root
                break
                
        if bin_dir:
            print(f"Poppler extracted to {bin_dir}")
            print("Adding to current session PATH...")
            os.environ["PATH"] += os.pathsep + bin_dir
            
            # Permanently add to User PATH using Windows Registry
            print("Attempting to add Poppler to your User PATH permanently...")
            try:
                import winreg
                with winreg.OpenKey(winreg.HKEY_CURRENT_USER, r'Environment', 0, winreg.KEY_ALL_ACCESS) as key:
                    current_path, _ = winreg.QueryValueEx(key, 'PATH')
                    if bin_dir not in current_path:
                        new_path = current_path + ';' + bin_dir
                        winreg.SetValueEx(key, 'PATH', 0, winreg.REG_EXPAND_SZ, new_path)
                        
                        # Notify the system of the environment change
                        import ctypes
                        HWND_BROADCAST = 0xFFFF
                        WM_SETTINGCHANGE = 0x001A
                        SMTO_ABORTIFHUNG = 0x0002
                        ctypes.windll.user32.SendMessageTimeoutW(
                            HWND_BROADCAST, WM_SETTINGCHANGE, 0, 'Environment', SMTO_ABORTIFHUNG, 5000, ctypes.byref(ctypes.c_long())
                        )
                        print("Successfully added Poppler to your system PATH.")
                        print("NOTE: You may need to restart your terminal/IDE for the PATH changes to take effect.")
                    else:
                        print("Poppler is already in your system PATH.")
            except Exception as e:
                print(f"Failed to add Poppler to permanent PATH: {e}")
                print(f"Please manually add this directory to your system PATH: {bin_dir}")
        else:
            print("Error: Could not find pdftoppm.exe in the extracted poppler files.")
    else:
        print(f"Unsupported OS: {os_name}. Please install poppler manually.")

def main():
    print("--- Python Packages ---")
    install_python_packages()
    
    print("\n--- System Dependencies ---")
    if check_poppler():
        print("Poppler is already installed and in PATH.")
    else:
        install_poppler()
        
    print("\nDependency check complete! You can now run the figure generation scripts.")

if __name__ == "__main__":
    main()
