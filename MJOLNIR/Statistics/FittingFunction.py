import matplotlib.pyplot as plt
import numpy as np


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
    
highlighter = '\\boldsymbol{'
ender= '}'
def executable(func):
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
        self.variableNames = [self.format[i]+self.parameterNames[i] for i in range(len(self.format))]
    @executable
    def setParameter(self,event,index):
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
        string = '$f(x, '
        string+=', '.join([x for x in self.variableNames])+') = '
        string+='{}'.format(self.variableNames[0])+'\\cdot \mathrm{exp}^{\left(-\\frac{\left(x-'+'{}'.format(self.variableNames[1])\
        +'\\right)^2}{2\\cdot '+'{}'.format(self.variableNames[2])+'^2}\\right)}+'+'{}'.format(self.variableNames[3])+'$'
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
        self.variableNames = [self.format[i]+self.parameterNames[i] for i in range(len(self.format))]
        self.currentFitParameter = 0
    @executable
    def setParameter(self,event,index):
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
        string = '$f(x, '
        string+=', '.join([x for x in self.variableNames])+') = '
        string+='{}'.format(self.variableNames[0])+'\\cdot\\frac{'+'{}'.format(self.variableNames[2])\
        +'^2}{\left(x-'+'{}'.format(self.variableNames[1])+'\\right)^2+'+'{}'.format(self.variableNames[2])+'^2}+'\
        +'{}'.format(self.variableNames[3])+'$'
        
        if not highlight is None:
            for i in highlight:
                self.variableNames[i] = self.variableNames[i][len(highlighter):-len(ender)]

        return string