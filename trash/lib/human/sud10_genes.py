import lib.genes
import lib.human.genes
import lib.expression

class GeneOrientatedPeaks(lib.human.genes.GeneOrientatedPeaks):
  def __init__(self, type):
    super(GeneOrientatedPeaks, self).__init__(type)
   
    # Add the extra annotation
   
    
