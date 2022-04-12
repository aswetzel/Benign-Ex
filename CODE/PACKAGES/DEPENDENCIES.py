import sys, os, glob, site, importlib, subprocess
from collections import defaultdict

if len(sys.argv) != 2:
    print "DEPENDENCIES.py -R <R CMD>"
    sys.exit(1)

try:
 try:
  PYTHON_VER = sys.version.split('\n')[0].split(' ')[0].strip()  # Get Python Version
  R = sys.argv[1].strip()
  R_VER = sys.argv[1].strip() + ' --version'  # Get R Version from AUTOMATION File

  # Check if Required Python Modules are Installed
  MODULES = ['progress', 'networkx', 'intervaltree']
  PNOTINSTALLED = defaultdict(list)
  PNTINSTALLIST = []
  for enu, modu in enumerate(MODULES):
    try:
        importlib.import_module(modu.strip())
    except ImportError:
        if modu.strip() == 'networkx':
            PNOTINSTALLED[modu.strip() + '==1.8.1'] = ['NA', enu]
            PNTINSTALLIST.append(modu.strip() + '==1.8.1')
        else:
            PNOTINSTALLED[modu.strip()] = ['NA', enu]
            PNTINSTALLIST.append(modu.strip())

  # Check if Required R Libraries are Installed
  PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
  DEPENDENCIES = glob.glob(PATH_TO_SCRIPT + '/*.tar.gz') + glob.glob(PATH_TO_SCRIPT + '/*.zip')
  COMMAND = R + ' CMD BATCH ' + PATH_TO_SCRIPT + '/PACKAGE_DEP.R ' + PATH_TO_SCRIPT + '/PACKAGE_DEP.Rout'
  try:
   if os.system(COMMAND) != 0: raise Exception
  except :
   raise KeyboardInterrupt
  ROUT = open(PATH_TO_SCRIPT + '/PACKAGE_DEP.Rout', 'r').readlines()
  PACKAGES = [ ro for ro, ut in enumerate(ROUT) if 'Version' in ut.strip() or 'proc.time' in ut.strip() ]
  INSTALLED = [ pckg.strip().split(' ')[0].strip() for pckg in ROUT[PACKAGES[0]:PACKAGES[1]][1:-1] ]
  LIBRARIES = [ lib.split("/")[-1].split("_")[0].strip() for lib in glob.glob(PATH_TO_SCRIPT + '/*.tar.gz') + glob.glob(PATH_TO_SCRIPT + '/*.zip') ]
  RNOTINSTALLED = list(set(LIBRARIES).difference(INSTALLED))

  # Install required Python modules via 'pip'
  if len(PNTINSTALLIST) > 0:
    PATH = os.path.dirname(os.path.realpath('__file__'))
    PIPINSTALL = ''
    try:
        PIPINSTALL = PNOTINSTALLED['pip'][0].strip()
    except IndexError:
        PIPINSTALL = 'YES'

    if PIPINSTALL != 'YES': # If pip not installed, install it
        COMMAND = 'python ' + PATH + '/get-pip.py'
        os.system(COMMAND)
        reload(site)
        import pip
    else:
        import pip
        PNTINSTL = PNTINSTALLIST
        for entinsl in PNTINSTL:
            print 'INSTALLING MODULE... ', entinsl.strip()
            try:
		subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "--disable-pip-version-check", entinsl.strip()])
            	print 'Module... ', entinsl.strip(), ' was installed successfully \n'
            except ImportError:
            	print 'Module... ', entinsl.strip(), " didn't install. Please try installing it manually."
		sys.exit()
        reload(site)

  # Install required R packages
  if len(RNOTINSTALLED) > 0: 
    for epckg in RNOTINSTALLED:
	for each in DEPENDENCIES: 
           if each.split("/")[-1].split("_")[0].strip() == epckg:
              print 'Library: ', epckg.strip(), " doesn't exist. Please install it manually."
    sys.exit()
  # Report if All Python Modules & R Packages are Installed
  else:
    print 'ALL THE DEPENDENCIES ARE UP TO DATE...'

 except IndexError:
  raise KeyboardInterrupt
except KeyboardInterrupt:
 print ''
 sys.exit(1)
