Optimizations
=============
This is a small set of notes explaining the methods used to optimize different parts of the 
software code. It is meant as an overview of the changes and thoughts that have gone into 
the code structure and changes made.



.. toctree::
   :maxdepth: 2

   OptimizationOfVoronoiTessellation
   
   

Timing function
_______________
In order to quantify whether or not a speed-up has been achieved, one need to time the 
method in question. This is mostly easily done through the use of the package *time* 
from Python, where one simply requests the current time in seconds before and after 
the method has been running. Subtracting the two then gives an estimate of the time 
consumption. The reason for the use of the word 'estimate' is that this is not the 
correct time used by the computer. The true time is the actual CPU time or even 
better, the number of CPU cycles needed to run the method. The reason the *time* module 
does not capture this is that if the program does not run as the only process and with 
a 100% usage of the core then a discrepancy between measured and actual time consumption 
is created. However, these technicalities are not taken into account when I have 
performed the profiling of the code as I in no case am able to perform a superb 
optimization with a background in pure physics.

The decorator used for profiling is given below:

.. code-block:: python
	:linenos:

	def my_timer_N(N=3):
	    if N<0:
		raise AttributeError('Number of runs need to be bigger or equal to 1 or equal to 0 for no timing, but {} given.'.format(N))
	    def my_timer(func):
		import time
		def newFunc(*args,**kwargs):
		    Time = []
		    if N ==0:
		        returnval = func(*args,**kwargs)
		    else:
		        for i in range(N):
		            startT = time.time()
		            returnval = func(*args,**kwargs)
		            stopT = time.time()
		            Time.append(stopT-startT)
		        if N>1:
		            print('Function "{}" took: {}s (±{}s)'.format(func.__name__,np.mean(Time),np.std(Time)/np.sqrt(N)))
		        else:
		            print('Function "{}" took: {}s'.format(func.__name__,Time[0]))
		    return returnval
		return newFunc
	    return my_timer

With this definition of the decorator, using it to time a function is straight forward, i.e. 

.. code-block:: python
	:linenos:
	
	
	def untimedFunction(*args,**kwargs):
		# some calculations
		return True
		
	@my_timer_N(N=5)
	def timedFunction(*args,**kwargs):
		# some calculations
		return True
	

Thus, the above code allows for timing the function *timedFunction* called 5 times, 
producing the output: Function "timedFunction" took: 1.1920928955078125e-06s 
(±1.3571023436315258e-06s), where the first time is the average of the N=5 runs, 
and the parenthesis denotes the uncertainty on the mean, :math:`std(X)/\sqrt{N}`.
	

	




   
   
   
