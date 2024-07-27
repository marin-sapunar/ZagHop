import numpy as np

class ElStates():
    def __init__(self, states, parent=None):
        self._parent = parent
        self._groups = None
        if isinstance(states, list):
            self._groups = [ElStates(group, parent=self) for group in states]
            states = {}
        self._sym = states.pop("symmetry", None)
        self._mult = states.pop("multiplicity", None)
        self._nstate = states.pop("nstate", 0)

    @property
    def sym(self):
        if self._sym is not None:
            return self._sym
        if self._parent is not None:
            return self._parent.sym
        return None

    @property
    def mult(self):
        if self._mult is not None:
            return self._mult
        if self._parent is not None:
            return self._parent.mult
        return None

    @property
    def nstate(self):
        if self._groups is not None:
            return np.sum([grp.nstate for grp in self._groups])
        return self._nstate

    def __iter__(self):
        if self._groups is not None:
            for group in self._groups:
                for i in group:
                    yield i
        else:
            yield self

    def __getitem__(self, key):
        if self._groups is not None:
            return self._groups[key]
        else:
            if key == 0:
                return self
            raise IndexError("list index out of range")
