import matplotlib.pyplot as plt
import numpy as np

from distutils.spawn import find_executable
if find_executable('latex'):
    USETEX = True
else:
    USETEX = False


class FittingFunction(object):
    def __init__(self,function):
        self.function = function
        self.__name__ = type(self).__name__
        self.executable = False
        self.fitError = False
    
    def __call__(self, x):
        return self.function(x,*self.parameters)

    @property
    def function(self):
        return self._function

    @function.getter
    def function(self):
        return self._function

    @function.setter
    def function(self,function):
        self._function = function
        self.parameterLength = self.function.__code__.co_argcount-2
        self.parameters = np.array([np.NaN]*self.parameterLength)
        self.parameterNames = np.array(self._function.__code__.co_varnames[2:])

    
    @property
    def parameters(self):
        return self._parameters

    @parameters.getter
    def parameters(self):
        return self._parameters

    
    @parameters.setter
    def parameters(self,parameters):
        self._parameters = parameters
        if np.sum(np.isnan(self._parameters))==0:
            self.executable=True
        else:
            self.executable=False
        
    

    def setParameter(self,event,index):
        raise NotImplementedError('The setParameter()-method is not yet implemented!')
    

highlighter = r'\left\{ '
ender= r'\right\} '
def executable(func): # pragma: no cover
    def newFunc(self,*args,**kwargs):
        returnval = func(self,*args,**kwargs)
        if np.sum(np.isnan(self.parameters))==0:
            self.executable=True
        else:
            self.executable=False
        return returnval
    return newFunc

class Gaussian(FittingFunction):
    def func(self,x,A,mu,sigma,B):
        return A * np.exp(-np.power(x-mu,2.0)/(2*sigma**2))+B
    
    def __init__(self):
        super(Gaussian,self).__init__(self.func)
        self.format=['','\\','\\','']
        #if USETEX == True:
        self.variableNames = [self.format[i]+self.parameterNames[i] for i in range(len(self.format))]
        #else:
        #    self.variableNames = [self.parameterNames[i] for i in range(len(self.format))]
    @executable
    def setParameter(self,event,index): # pragma: no cover
        if index == 0:
            if not np.isnan(self.parameters[3]):
                self.parameters[0] = event.ydata-self.parameters[3]
            else:
                self.parameters[0] = event.ydata
            self.parameters[1] = event.xdata
            return 2
        elif index == 1:
            self.parameters[1] = event.xdata
            return 2
        elif index == 2:
            self.parameters[2] = np.abs(self.parameters[1]-event.xdata)/(np.sqrt(2*np.log(2)))
            return 3
        elif index == 3:
            if np.isnan(self.parameters[3]):
                self.parameters[0]-= event.ydata
            else:
                self.parameters[0]-= event.ydata-self.parameters[3]
            self.parameters[3] = event.ydata
            return 0


    def latex(self,highlight=None):
        if not highlight is None:
            if highlight==0:
                highlight = [0,1]
            highlight = np.array(highlight).reshape(-1)
            for i in highlight:
                self.variableNames[i] = highlighter+self.variableNames[i]+ender
        #if USETEX==True:
        string = '$f(x, '
        string+=', '.join([x for x in self.variableNames])+') = '
        string+='{}'.format(self.variableNames[0])+'\\cdot \mathrm{exp}^{\left(-\\frac{\left(x-'+'{}'.format(self.variableNames[1])\
        +'\\right)^2}{2\\cdot '+'{}'.format(self.variableNames[2])+'^2}\\right)}+'+'{}'.format(self.variableNames[3])+'$'
        #else:
        #    string = '$f(x, '
        #    string+=', '.join([x for x in self.variableNames])+') = '
        #    string+='{}'.format(self.variableNames[0])+'*exp(-[x-'+'{}'.format(self.variableNames[1])\
        #    +']^2/[2*'+'{}'.format(self.variableNames[2])+'^2])+'+'{}'.format(self.variableNames[3])+'$'
        if not highlight is None:
            for i in highlight:
                self.variableNames[i] = self.variableNames[i][len(highlighter):-len(ender)]
        return string
    

class Lorentz(FittingFunction):
    def func(self,x,A,x_0,gamma,B):
        return A * np.power(gamma,2.0)/(np.power(x-x_0,2.0)+np.power(gamma,2.0))+B
    
    def __init__(self):
        super(Lorentz,self).__init__(self.func)
        self.format = ['','','\\','']    
        #if USETEX == True:
        self.variableNames = [self.format[i]+self.parameterNames[i] for i in range(len(self.format))]
        #else:
        #    self.variableNames = [self.parameterNames[i] for i in range(len(self.format))]
        self.currentFitParameter = 0
    @executable
    def setParameter(self,event,index): # pragma: no cover
        if index == 0:
            self.parameters[0] = event.ydata
            self.parameters[1] = event.xdata
            return 2
        elif index == 1:
            self.parameters[1] = event.xdata
            return 2
        elif index == 2:
            self.parameters[2] = np.abs(self.parameters[1]-event.xdata)
            return 3
        elif index == 3:
            if np.isnan(self.parameters[3]):
                self.parameters[0]-= event.ydata
            else:
                self.parameters[0]-= event.ydata-self.parameters[3]
            self.parameters[3] = event.ydata
            return 0     

    def latex(self,highlight=None):
        if not highlight is None:
            if highlight==0:
                highlight = [0,1]
            highlight = np.array(highlight).reshape(-1)
            for i in highlight:
                self.variableNames[i] = highlighter+self.variableNames[i]+ender
        #if USETEX == True:
        string = '$f(x, '
        string+=', '.join([x for x in self.variableNames])+') = '
        string+='{}'.format(self.variableNames[0])+'\\cdot\\frac{'+'{}'.format(self.variableNames[2])\
        +'^2}{\left(x-'+'{}'.format(self.variableNames[1])+'\\right)^2+'+'{}'.format(self.variableNames[2])+'^2}+'\
        +'{}'.format(self.variableNames[3])+'$'
        #else:
        #    string = 'f(x, '
        #    string+=', '.join([x for x in self.variableNames])+') = '
        #    string+='{}'.format(self.variableNames[0])+'*('+'{}'.format(self.variableNames[2])\
        #    +'^2)/([x-'+'{}'.format(self.variableNames[1])+']^2+'+'{}'.format(self.variableNames[2])+'^2)+'\
        #    +'{}'.format(self.variableNames[3])
        if not highlight is None:
            for i in highlight:
                self.variableNames[i] = self.variableNames[i][len(highlighter):-len(ender)]

        return string


class InclineLinear(FittingFunction):
    def func(self,x,a,b,c):
        d = (c-b)/a
        y = np.zeros_like(x)
        y[x<d] = a*x[x<d]+b
        y[x>=d] = c
        return y
    
    def __init__(self):
        super(InclineLinear,self).__init__(self.func)
        self.format=['','','']
        #if USETEX == True:
        self.variableNames = [self.format[i]+self.parameterNames[i] for i in range(len(self.format))]
        #else:
        #    self.variableNames = [self.parameterNames[i] for i in range(len(self.format))]
    @executable
    def setParameter(self,event,index): # pragma: no cover
        if index == 0:
            self.parameters[1] = event.ydata
            return 1
        elif index == 1:
            if not np.isnan(self.parameters[1]):
                self.parameters[0] = (event.ydata-self.parameters[1])/event.xdata
            else:
                self.parameters[0] = event.ydata/event.xdata
            return 2
        elif index == 2:
            self.parameters[2] = event.ydata#np.abs(self.parameters[1]-event.xdata)/(np.sqrt(2*np.log(2)))
            return 0
        


    def latex(self,highlight=None):
        if not highlight is None:
            if highlight==0:
                highlight = 1
            elif highlight == 1:
                highlight = 0
            highlight = np.array(highlight).reshape(-1)
            for i in highlight:
                self.variableNames[i] = highlighter+self.variableNames[i]+ender
        #if USETEX==True:
        string = r'$f(x, '
        string+=r', '.join([x for x in self.variableNames])+') = '
        string+=r'\left\{ \stackrel{'+str(self.variableNames[0])+r'x+'+str(self.variableNames[1])+r' \quad x<\frac{'+str(self.variableNames[2])\
        +r'-'+str(self.variableNames[1])+r'}{'+str(self.variableNames[0])+r'}}{ '+str(self.variableNames[2])+r' \quad x\geq\frac{'+str(self.variableNames[2])\
        +r'-'+str(self.variableNames[1])+r'}{'+str(self.variableNames[0])+r'}} \right. $'
        #else:
        #    string = 'f(x, '
        #    string+=', '.join([x for x in self.variableNames])+') = '\
        #    +'if x<('+str(self.variableNames[2])+'-'+str(self.variableNames[1])+')/'+str(self.variableNames[0])+': '+str(self.variableNames[0])+'x+'\
        #    +str(self.variableNames[1])+' else '+str(self.variableNames[2])
        if not highlight is None:
            for i in highlight:
                self.variableNames[i] = self.variableNames[i][len(highlighter):-len(ender)]
        return string










def test_FittingFunction_Call():
    for USETEX in [True,False]:
        for f in FittingFunction.__subclasses__():
            fun = f()
            parameters = np.random.rand(fun.parameterLength)
            fun.parameters = parameters
            assert(fun.executable==True)

            x = np.linspace(-5,5,201)
            y = fun(x)

            text = fun.latex(highlight=0)
