import sys
import math
import random
import matplotlib.pyplot as plt
import numpy.fft as fft
import numpy as np
import scipy.integrate as integrate


#библиотеки




def values():
    filename = sys.argv[1]
    i = 0
    f = open(filename)
    x = []  #временная полоса
    y = [] #мощность
    for elem in f:
        x_i=elem.split(";")[0]
        y_i= elem.split(";")[1]
        if i == 0:  # скип первой строчки
            i = 1
            continue
        y.append(int(y_i))
        x.append(int(x_i))
    return [x,y]


def Lin_Aprox(xy):  # МНК https://prog-cpp.ru/mnk/
    n=len(xy[0])
    sumx = sum(xy[0])
    sumy = sum(xy[1])
    sumxy=0
    sumx2=0
    for i in range (n):
        sumxy+=xy[0][i]*xy[1][i]
        sumx2+=pow(xy[0][i],2)
    a=(n*sumxy-sumx*sumy)/(n*sumx2-sumx*sumx)
    b=(sumy-a*sumx)/n
    return (a,b)


def Noise(xy,ab): #реализация шума (эпсилент от t)
    x_n=xy[0] #иксовая координата естественно без шума
    y_n=[]
    for i in range (max(xy[0])):
        y_n.append(xy[1][i]-ab[0])

    return [x_n, y_n]


def hist_of_noise(xy_noise): #получить значения инервалов и соотвествующие количество точек в этих интервалах
                                #для шума, шаг условно равен 5
    Min=min(xy_noise[1])
    Max=max(xy_noise[1])

    l=round((Max-Min)/5)

    y_hist = [i * 0 for i in range(0, l)] #количество точек в отрезках ,которым принадлежат середины


    h=5
    x_hist =  [ Min+h*i-h/2 for i in range(0,l) ]    #координаты середин отрезков
    for elem in xy_noise[1]:
        b=int((elem-Min)/h)
        if b==16:
            b=b-1
        y_hist[b] += 1

    return [x_hist,y_hist]



def GraphOfTimeRow(xy):  #  график временного ряда
    plt.plot(xy[0], xy[1],'black')





def Spectr(xy):  #график амплитудного спектра
#это формула из лекции по преобразователям информации, длину делют на 2 для её стабилизации амплитуы
    ft = (abs(fft.fft(xy[1], len(xy[1])) / (len(xy[1]) / 2)))  # Используем БПФ
    x=[(i+1)/ len(xy[1]) for i in range( len(xy[1]))]
    plt.plot(x, ft)
    plt.title("Амплитудный спектр сигнала")
    plt.xlabel("Частота, 1/отн.ед.")
    plt.ylabel("СП амплитуды, Вт*отн.ед.")
    plt.axis([0, 1, -0, 2.5])
    plt.grid()
    plt.show()




def PlotAprox(xy,ab):  # Изображения функции линейной аппроксимации сигнала
    x_line = xy[0]
    y_line = [elem*ab[0]+ab[1] for elem in x_line]
    plt.plot(xy[0], xy[1], 'black')
    plt.plot(x_line,y_line,'r')
    plt.axis([0, 1024, 360, 480])
    plt.title("Линейная аппроксимация сигнала")
    plt.grid()
    plt.xlabel("Время, отн.ед.")
    plt.ylabel("Мощность, Вт")
    plt.show()
    return [x_line,y_line]


def plotNoise(n_xy):  # График шумовой компоненты
    plt.plot(n_xy[0], n_xy[1])
    plt.title("График шумовой зависимости компоненты E(t)")
    plt.grid()
    plt.xlabel("Время, отн.ед.")
    plt.ylabel("Мощность, Вт")
    plt.show()

def plotHist(hist_n_xy):  #  гистограмма шума
    hist_n_xy1=[elem - 416.9 for elem in hist_n_xy[0]]      #число 416.9 середина области абцисс , и вычитать надо
    # так как мы рассматриваем шум , очень сложно описать почему надо вычитать
    # но чутсвуется что так надо, инае формула создающая
    #апроксимирующую функцию будет возвращать нули
    plt.plot(hist_n_xy1,hist_n_xy[1],'ob')
    x_appr =[max(hist_n_xy1)+min(hist_n_xy1) - min(hist_n_xy1)*i/100  for i in range(-100,101)]


    y_appr = [probability_density_function(elem) for elem in x_appr]
   # print(y_appr)
    #насколько я понял это построение апроксимирующей гистограмму функции
    plt.plot(x_appr,y_appr, 'r')
    plt.axis([-40,40,0,140])
    plt.grid()
    plt.xlabel("Мощность, Вт")
    plt.ylabel("Количество точек, шт.")
    plt.title("Гистограмма шумовой компоненты E(t)")

    plt.show()

def probability_density_function(x):  # вспомогательная функция из текста лабы
    n = 5
    y_0 = 0
    x_0 = 0
    sigma = 22
    A = 120
    ans = y_0 + A * math.exp(-(((x - x_0) / (sigma)) ** (2 * n)))
    return ans



def drawPredict(xy_apr,ab):  # Изображение области предсказания
    #https: // coderlessons.com / tutorials / python - technologies / uchitsia - stsipi / scipy - integratsiia
    #как использовать функцию отсюда
    x_appr = xy_apr[0]
    y_appr = xy_apr[1]
    for i in range(1, 21):
        x_appr.append(1000 + i)
        y_appr.append(0.0014991434991434992 * ab[0]*x_appr[1000 + i - 1] + 424.1366786786787)
    mean_val = integrate.quad(help_1, -np.inf, np.inf)#интегрирование
    print(help_1)
    std = math.sqrt(integrate.quad(help_2, -np.inf, np.inf)[0]) - mean_val[0] ** 2
    print("велечина отклонения за счёт шумовой составляющей =",std,"ВТ")
    print("точки прогноза будут лежать в диапозоне= [{0} , {1}] ВТ".format(-std,std))
    ANS = integrate.quad(help_3,-1*std, std)
    print(" с вероютностью=",ANS,"%")

    plt.plot(x_appr, y_appr, 'r')
    line1 = [elem - std for elem in y_appr]
    line2 = [elem + std for elem in y_appr]
    plt.plot(x_appr,line1 , 'b')
    plt.plot(x_appr, line2, 'b')
    plt.axis([1000, 1024, 400, 450])
    plt.title("Область, в которой с вероятностью 60% будут лежать результаты")
    plt.grid()
    plt.xlabel("Время, отн.ед.")
    plt.ylabel("Мощность, Вт")
    plt.show()



help_1= lambda x:0 + 120 * math.exp(-((x - 0) / (22)) ** (2 * 5)) * x
help_2= lambda x:0 + 120 * math.exp(-((x - 0) / (22)) ** (2 * 5)) * x*x/5023.132
help_3= lambda x:0 + 120 * math.exp(-((x - 0) / (22)) ** (2 * 5))/5023.132







if __name__ == "__main__":

    xy = values()
    ab=Lin_Aprox(xy)

    n_xy = Noise(xy, ab)
    hist_n_xy = hist_of_noise(n_xy)


    GraphOfTimeRow(xy)
    Spectr(xy)
    xy_apr = PlotAprox(xy,ab)
    plotNoise(n_xy)
    plotHist(hist_n_xy)
    drawPredict(xy_apr,ab)