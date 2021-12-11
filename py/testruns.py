import numpy as np
import perple_x as run
import parameters as p

# stars = ["HIP 98355", "sun"]
run.start(name='marstest',
          M_p=0.1*p.M_E, stars='sun', core_efficiency=0.8,  # planet parameters
          n=1000, tol=0.1,  # interior structure numerical parameters
          overwrite=True, verbose=False  # run settings
          )


