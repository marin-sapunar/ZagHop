""" Base class for interfacing with electronic structure codes. """
import sys


class QMInterface():
    """ 
    Generic interface for electronic structure calculations.

    The main methods of this class are `update`, `run`, and `read`. The
    `update` method is used to set up the calculation, `run` is used to execute
    the calculation, and `read` is used to extract the requested results from
    the output files.

    Derived classes should implement the following methods:
    - `check_template`: Read the input file and save any information that will
         be required by the other methods. Additionally, it should define the
         read_funcs dictionary.
    - `run`: Execute the QM calculation.
    - `update`: Update all inputs based on the provided information on the 
                system and requested values.
    - `read_*` : Methods to be called to read values from the QM output.

    Methods starting with `read_` are called by the `read` method and should
    store the calculated values in the `results` dictionary.
        
    Attributes:
        read_funcs (dict) : Key, value pairs of values which can be read by
                            the interface and the corresponding read functions.
        options (dict): Options and paths relevant for the QM calculation.
        system (dict): Information on the system.
        request (list): Values which should be calculated/read.
        results (dict): Results read from the QM output.
        calc_done (bool): True after a calculation has been completed.
    """
    def __init__(self,
                 template="qm.inp",
                 log_file="qm.log",
                 err_file="qm.err",
                 **kwargs):
        """
        Initialize the QM interface object.
        
        Parameters:
            - template (str): Path to the template file. Default is "qm.inp".
            - log_file (str): Path to the log file. Default is "qm.log".
            - err_file (str): Path to the error file. Default is "qm.err".
            - **kwargs: Additional keyword arguments.
        """
        self.calc_done = False
        self.read_funcs = {}
        self.system = {}
        self.results = {}
        self.request = {}
        self.options = {
            "template": template,
            "log_file": log_file,
            "err_file": err_file
        }
        self.options.update(kwargs)

    @classmethod
    def generate_inputs(cls, ipath, system, options):
        """ Generate input templates to be used by the interface. """
        raise NotImplementedError("Need to call specific interface.")

    def check_template(self):
        """ Read interface options from input files. """
        raise NotImplementedError("Need to call specific interface.")

    def update(self, system, request):
        """ 
        Update input files to prepare a calculation. 
        
        Parameters:
            - system: Information on the system.
            - request: Values to be calculated/read.
        """
        raise NotImplementedError("Need to call specific interface.")

    def run(self):
        """ Run the calculation and check success. """
        raise NotImplementedError("Need to call specific interface.")

    def read(self):
        """ Read all requested values.
        
        Returns:
            The data dictionary containing the calculated values.
        """
        for req in self.request:
            if req not in self.read_funcs:
                print(f"Error in QM interface. {req} not available.")
                sys.exit(1)
            self.read_funcs[req]()
        return self.results

    def __getattr__(self, atr):
        """ Convenience for easier access to options/system variables."""
        try:
            return self.options[atr]
        except KeyError:
            return self.system[atr]
