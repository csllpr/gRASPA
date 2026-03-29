import os
import shutil
import subprocess

homedir = os.getcwd()
repodir = os.path.dirname(homedir)
basics = ['CO2-MFI', 'Methane-TMMC', 'Bae-Mixture', 'NU2000-pX-LinkerRotations', 'Tail-Correction']
ref_calc = ['Reference_NIST_SPCE/Box-1/', 'Reference_NIST_SPCE/Box-2/', 'Reference_NIST_SPCE/Box-3/', 'Reference_NIST_SPCE/Box-4/']
sims = []
sims.extend(basics)
sims.extend(ref_calc)

binary = os.path.join(repodir, "src_clean", "nvc_main.x")
cleanup_dirs = ["AllData", "FirstBead", "Lambda", "Movies", "Restart", "TMMC"]

for direct in sims:
  os.chdir(direct)
  with open("output.txt", "w") as output_file:
    subprocess.run([binary], stdout=output_file, stderr=subprocess.STDOUT, check=True)
  for dirname in cleanup_dirs:
    shutil.rmtree(dirname, ignore_errors=True)
  os.chdir(homedir)
  print(f"Simulation {direct} has finished.\n")
