import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore

def version():
    return "ykktmm 3 [24/06/24]"

DEFAULT = object()
rho0 = 1.21
c0 = 343

#Numa proxima versao apenas
class templates:
    def tubo(self):
        return
    def barrier(self):
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
    def __init__(self,area,tam=0,zc=DEFAULT,kc=DEFAULT,assoc=">"): #Tem muita coisa errada aqui
        self.area = area
        self.tam = tam
        self.zc = zc
        self.kc = kc
        self.assoc = assoc

    #Atualiza zc e kc para campo livre caso a variavel nao tenha sido inicializada
    def updateProperties(self):
        if self.zc is DEFAULT:
            self.zc = rho0*c0*np.ones(self.freq.shape)
        if self.kc is DEFAULT:
            self.kc = 2*np.pi*self.freq/c0

        
    def setFreqVec(self,f):
        self.freq = f
        self.updateProperties()

    #getters e setters servem para facilitar a expansão do código
    def getMatrix(self):
        one_vector = np.ones(self.freq.shape) #Força a tudo ser do mesmo tamanho na direcao da frequencia
        zero_vector = np.zeros(self.freq.shape) #Otimizando a performance

        if (self.assoc == "series") or self.assoc == ">":
            if self.tam == 0:
                print("Utilize o argumento tam para definir o comprimento da camada")
                print("A camada será ignorada(Matriz identidade)")
            arg = self.kc*self.tam
            matrix = np.array([[np.cos(arg)   ,   1j*(self.zc/self.area)*np.sin(arg)],
                            [1j*(self.area/self.zc)*np.sin(arg) ,   np.cos(arg)]])

        elif (self.assoc == "parallel") or self.assoc == "//": #Eventualmente renomear
            matrix = np.array([[ one_vector  ,    zero_vector            ],
                            [ one_vector*(self.area/self.zc)  , one_vector]])
        
        elif (self.assoc == "barrier") or self.assoc == "||":
            matrix = np.array([[ one_vector  ,    one_vector*(self.zc/self.area)],
                            [zero_vector ,       one_vector               ]])
            
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
                A forma que eu montei o código resultou em uma estrutura de dados meio esquisita
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
        
        else:
            return


    def solveZs(self,zl=DEFAULT):
        tg = self.solveTg()
        #Paredes Rígidas
        if zl is DEFAULT:
            zs = self.layerlist[0].area * (tg[0,0]/tg[1,0])

        #Impedância de radiação conhecida
        else:
            zs = tg[0,1] + (tg[0,0]*zl)/self.layerlist[-1].area
            zs /= tg[1,1] + (tg[1,0]*zl)/self.layerlist[-1].area
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
        return 1 - (np.abs((zs-rho0*c0)/(zs+rho0*c0)))**2
    
    def dB(val, ref=1):
        return 20*np.log10(abs(val)/ref)

    def delany_basley(f,sigma):
        #Ref: Cox pg 141
        X = (rho0*f)/sigma
        zc = rho0*c0*(1 + 0.0571*(X**-0.754) - 1j*0.087*(X**-0.732))
        kc = ((2*np.pi*f)/c0)*(1 + 0.0978*(X**-0.700) + -0.189j*(X**-0.595))
        return zc,kc
    
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