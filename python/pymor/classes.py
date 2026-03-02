from dune.istl._istl import BlockVector as DuneVector

from dune.darcyflow import _darcyflow as df

from pymor.vectorarrays.list import CopyOnWriteVector

class WrappedDuneVector(CopyOnWriteVector):
  """
  Wrapper class for DUNE-internal Vector with pyMOR IF
  """

  def __init__(self, vector):
    assert isinstance(vector, DuneVector)
    self._impl = vector

  
