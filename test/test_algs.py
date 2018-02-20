import numpy as np
import smithWater  

def testSmithWaterman():
    """
    check if pos pair scores are greater than neg pair scores 
    """
    assert smithWater.getAverageScores("MATIO", -7, -3)[0] > smithWater.getAverageScores("MATIO", -7, -3)[1]
    assert smithWater.getAverageScores("BLOSUM50", -7, -3)[0] > smithWater.getAverageScores("BLOSUM50", -7, -3)[1]
    

    
