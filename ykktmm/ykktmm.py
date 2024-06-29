"""
ykktmm 3[edição Lista 2 de controle de ruido] (26/06/24)
Alterações:
    * air_props agora tem varios props
    * JCA()
"""

import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore

def version():
    return "ykktmm 3 [24/06/24]"

DEFAULT = object()

class air_props:
    #densidade do ar
    _rho0 = 1.21
    _c0 = 343
    #numero de prandtl
    _prandtl = 0.71
    _eta = 1.82E-5
    #cpcv
    _gamma = 1.4
    #pressao atmosferica
    _p0 = 101320

    @staticmethod
    def rho0(new_rho0 = DEFAULT):
        if new_rho0 is DEFAULT:
            return air_props._rho0
        else:
            air_props._rho0 = new_rho0
    
    @staticmethod
    def c0(new_c0 = DEFAULT):
        if new_c0 is DEFAULT:
            return air_props._c0
        else:
            air_props._c0 = new_c0
    
    @staticmethod
    def prandtl(new_prandtl = DEFAULT):
        if new_prandtl is DEFAULT:
            return air_props._prandtl
        else:
            air_props._prandtl = new_prandtl
            
    @staticmethod
    def eta(new_eta = DEFAULT):
        if new_eta is DEFAULT:
            return air_props._eta
        else:
            air_props._eta = new_eta
    
    @staticmethod
    def gamma(new_gamma = DEFAULT):
        if new_gamma is DEFAULT:
            return air_props._gamma
        else:
            air_props._gamma = new_gamma
    
    @staticmethod
    def p0(new_p0 = DEFAULT):
        if new_p0 is DEFAULT:
            return air_props._p0
        else:
            air_props._p0 = new_p0
    
    #retorna todas as propriedades definidas
    @staticmethod
    def current_props():
        return

"""
Implementação da classe camada
"""
class layer:

    """
    Attributes:
    area: Área em m^2
    tam: Tamanho L da camada
    zc: Impedância do Meio (a ser calculada por fora)
    kc: Número de onda do Meio (a ser calculado por fora)
    assoc: Tipo de associação (por padrão é em série)
    """
    def __init__(self,area=1,tam=0,zc=DEFAULT,kc=DEFAULT,assoc=">"):
        self.area = area
        self.tam = tam
        self.zc = zc
        self.kc = kc
        self.assoc = assoc

    #Atualiza zc e kc para campo livre caso a variavel nao tenha sido inicializada
    def updateProperties(self):
        if self.zc is DEFAULT:
            self.zc = air_props.rho0()*air_props.c0()*np.ones(self.freq.shape)
        if self.kc is DEFAULT:
            self.kc = 2*np.pi*self.freq/air_props.c0()

        
    def setFreqVec(self,f):
        self.freq = f
        self.updateProperties()


    def getMatrix(self):
        one_vector = np.ones(self.freq.shape) #Força a tudo ser do mesmo tamanho na direcao da frequencia
        zero_vector = np.zeros(self.freq.shape) #Otimizando a performance

        if (self.assoc == "series") or self.assoc == ">":
            if self.tam == 0:
                print("Utilize o argumento tam para definir o comprimento da camada")
                print("A camada será ignorada(Matriz identidade)")
            arg = self.kc*self.tam
            matrix = np.array([
                [np.cos(arg)   ,   1j*(self.zc/self.area)*np.sin(arg)],
                [1j*(self.area/self.zc)*np.sin(arg) ,   np.cos(arg)]
                ])

        elif (self.assoc == "parallel") or self.assoc == "//": #Eventualmente renomear
            matrix = np.array([
                [ one_vector  ,    zero_vector            ],
                [ one_vector*(self.area/self.zc)  , one_vector]
                ])
        
        elif (self.assoc == "barrier") or self.assoc == "||":
            matrix = np.array([
                [one_vector  ,    one_vector*(self.zc/self.area)],
                [zero_vector ,       one_vector               ]
                ])
            
        return matrix



class material:
    """
    Attributes:
    layerlist: Lista de todas as camadas e a forma de associação entre elas
    boundary: Condições de contorno para p0,q0.
    """
    def __init__(self,layers=[],boundary = DEFAULT,freq = []):
        self.layerlist = layers #muda de nome pq sim
        self.boundary = boundary
        self.freq = freq

        for layer in self.layerlist:
            layer.setFreqVec(freq)


    #Define as condições de contorno
    def setBound(self,p0,q0):
        # Caso seja independente da frequência, criar um puta vetor independente da frequência
        if not (type(p0) == np.ndarray or type(p0) == list):
            p0 = p0*np.ones(self.freq.shape)
        if not (type(q0) == np.ndarray or type(q0) == list):
            q0 = q0*np.ones(self.freq.shape)
        self.boundary = np.array([p0,q0])


    #Retorna a matriz de transferencia resultante
    def solveTg(self):
        tg = self.layerlist[0].getMatrix()
        n_tg = tg
        for i in range(1,len(self.layerlist)):
            tn = self.layerlist[i].getMatrix()
            #matmul manualmente para permitir a multiplicação termo a termo dos valores dependendo da frequencia
            n_tg[0,0] = tg[0,0]*tn[0,0] + tg[0,1]*tn[1,0]
            n_tg[0,1] = tg[0,0]*tn[0,1] + tg[0,1]*tn[1,1]
            n_tg[1,0] = tg[1,0]*tn[0,0] + tg[1,1]*tn[1,0]
            n_tg[1,1] = tg[1,0]*tn[0,1] + tg[1,1]*tn[1,1]
            tg = n_tg
        
        return tg


    #Retorna um vetor com [PL,QL], acho que ta quebrado #soacho
    def solveEq(self):
        if self.isBoundarySet():
            tg = self.solveTg()
            result = np.zeros([tg[0,0].shape[0],2],dtype=complex)
            for i in range(tg[0,0].shape[0]):
                """
                A forma que eu montei o código resultou em uma estrutura meio esquisita
                Então eu tenho que manualmente remontar a matriz para cada frequência
                """
                tn = np.array([
                    [tg[0,0][i], tg[0,1][i]],
                    [tg[1,0][i], tg[1,1][i]]
                               ])
                bn = np.array([self.boundary[0][i], self.boundary[1][i]])

                #Caso haja singularidade, pular o solve()
                # if np.linalg.det(tn) == 0:
                #     print("há uma matriz singular, ignorando ela")
                #     continue

                result[i] = np.linalg.solve(tn,bn)
        
            return np.transpose(result)
        
        #Se as condições de contorno não foram declaradas
        else:
            return


    def solveZs(self,zl=DEFAULT):
        tg = self.solveTg()
        #Paredes Rígidas
        if zl is DEFAULT:
            # print("zl eh default")    
            zs = self.layerlist[0].area * (tg[0,0]/tg[1,0])
        #Impedância de radiação conhecida
        else:
            # print("zl necas default")
            zs = tg[0,1] + (tg[0,0]*zl)/self.layerlist[-1].area
            zs /= (tg[1,1] + (tg[1,0]*zl)/self.layerlist[-1].area)
            zs *= self.layerlist[0].area
        
        return zs

    # Checa se as condições de contorno foram determinadas
    def isBoundarySet(self):
        if self.boundary is DEFAULT:
            print("Condições de contorno não específicadas")
            print("utilize setBound() ou as defina na inicialização do material")
            return False
        return True




class utils():
    def absort(zs):
        return 1 - (np.abs((zs-air_props.rho0()*air_props.c0())/(zs+air_props.rho0()*air_props.c0())))**2
    

    def dB(val, ref=1):
        #apenas realiza uma divisão termo a termo quando necessário
        if ref == 1:
            return 20*np.log10(np.abs(val))
        
        return 20*np.log10(np.abs(val)/ref)
        

    def delany_basley(f,sigma):
        #Ref: Cox pg 141
        X = (air_props.rho0()*f)/sigma
        zc = air_props.rho0()*air_props.c0()*(1 + 0.0571*(X**-0.754) - 1j*0.087*(X**-0.732))
        kc = ((2*np.pi*f)/air_props.c0())*(1 + 0.0978*(X**-0.700) + -0.189j*(X**-0.595))
        return zc,kc


    def JCA(f, sigma,tortuosidade, phi,comp_viscoso,comp_termico):
        omega = 2*np.pi*f
        eta = air_props.eta()
        gamma = air_props.gamma()
        Pr = air_props._prandtl()
        p0 = air_props.p0()
        rho0 = air_props.rho0()


        rhoef = np.sqrt( 1 + (4j*omega*rho0*eta*(tortuosidade**2))/((sigma*phi*comp_viscoso)**2) )
        rhoef = 1 + ((phi*sigma)/(1j*omega*rho0*tortuosidade)) * rhoef
        rhoef *= rho0*tortuosidade

        kef = np.sqrt(1 + 1j*(omega*Pr*rho0*(comp_termico**2))/(16*eta))
        kef = 1 - ( 8j*eta / (omega*Pr*(comp_termico**2)*rho0) )*kef
        kef = gamma - (gamma-1)/kef
        kef = gamma*p0/kef

        zc = np.sqrt(rhoef*kef)
        kc = omega*np.sqrt(rhoef/kef)
        return zc, kc
    
    def plot_cuxticks(minorticks=False,bandstyle="oct",showaudible=False,axis=plt): #Não me orgulho dessa função
        fr_minmax = axis.xlim()

        if bandstyle == 'oct':
            freq = [20,31.5,63,125,250,500,1000,2000,4000,8000,16000,20000]
            freq_txt = ["20","31.5","63","125","250","500","1k","2k","4k","8k","16k","20k"]


        if not showaudible:
            freq = freq[1:len(freq)-1]
            freq_txt = freq_txt[1:len(freq_txt)-1]

        ret_val = []
        ret_txt = []
        for i in range(len(freq)):
            if freq[i] < fr_minmax[0] or freq[i] > fr_minmax[1]:
                continue
            ret_val.append(freq[i])
            ret_txt.append(freq_txt[i])
        
        if not minorticks:
            axis.minorticks_off()

        #Testando isso aqui ainda
        #Inferno o animal que inventou a classe axis no matplotlib
        if axis is plt:
            axis.xticks(ret_val,ret_txt)
        else:
            axis.set_xticks(ret_val,ret_txt)
        return[ret_val,ret_txt]