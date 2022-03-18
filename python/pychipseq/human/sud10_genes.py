import pychipseq.genes
import pychipseq.human.genes
import pychipseq.expression

class GeneOrientatedPeaks(pychipseq.human.genes.GeneOrientatedPeaks):
  def __init__(self, type):
    super().__init__(type)
   
    # Add the extra annotation
   
    
