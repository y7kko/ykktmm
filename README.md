# ykktmm

O módulo tem objetivo de abstrair a maior parte da matemática envolvida no método das matrizes de transferência. Tornando a prototipagem de materiais um joguinho de lego.

## Sobre
#### O módulo foi feito para computar
* A impedância de superfície com incidência normal de um material multicamadas
* A matriz transferência resultante das camadas
* A pressão e velocidade de volume no final de um duto (NÃO TERMINADO)

#### Aprendendo a usar
Eu preparei um guia curto em Jupyter com exemplos. Ele se encontra [AQUI](https://github.com/y7kko/ykktmm/blob/main/Guia_Rapido.ipynb).

## Instalação
Para instalar basta executar o comando abaixo no seu console:

```
pip install git+https://github.com/y7kko/ykktmm.git
```

No Colab é a mesma coisa só que com um "!" antes:
```
!pip install git+https://github.com/y7kko/ykktmm.git
```

## Problemas(e incertezas)

* `solveEq()` apresenta valores esquisitos e incertos