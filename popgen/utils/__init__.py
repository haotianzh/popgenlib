import os
import popgen
import jpype
import jpype.imports
from jpype.types import *
from .treeutils import *
from .utils import *
package_dirname = os.path.dirname(popgen.__file__)
jpype.startJVM(classpath=[os.path.join(package_dirname, 'libs/*')])
from .javautils import *
from java.lang import System
from java.io import PrintStream, File
original = System.out
# System.setOut(PrintStream(File("NUL")))
