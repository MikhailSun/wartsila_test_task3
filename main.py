import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import makecharts as mc
from copy import copy

#Исходные данные (все единицы в СИ, если не указано иное)
#расчет квазинестационарный

#задаем напорную ветку характеристики компрессора. Компрессор всегда работает на постоянных оборотах.
compr_PR=[40,35,32,25,20,15,5,]
compr_G=[0.,1.,2.,3.,4.,4.5,5.]
compressor_curve_for_plot=interp1d(compr_G, compr_PR,bounds_error=False,kind='quadratic',fill_value='extrapolate')
compressor_curve_for_calculation=interp1d(compr_PR, compr_G,bounds_error=False,kind='quadratic',fill_value='extrapolate')
x_plot=np.linspace(0.,6.,50)
y_plot=[compressor_curve_for_plot(x) for x in x_plot]
fig1=mc.Chart(points_for_plot=[{'x':x_plot,'y':y_plot}],xlabel='G, kg/s',ylabel='Pressure Ratio', title="Compressor map", dpi=150,figure_size=(5,5))

eff_compr=0.7 #изоэнтропический КПД компрессора (для простоты примем его постоянным, но вообще это должна быть функция, как часть характеристики компрессора)
Cp=1050 #изобарная теплоемкость воздуха, Дж/кг/К (для упрощения принимаем постоянным, но вообще по-хорошему это должна быть функция от температуры)
k=1.4 #коэффициент адиабаты (для упрощения принимаем постоянным, но вообще по-хорошему это должна быть функция от температуры)
R=282. #Дж/кг/К, газовая постоянная для воздуха
P_atm=101325 #давление на выходе из всех клапанов (атмосферное)
T_compr_in=288.15 #температура воздуха на входе в компрессор
P_safety_valve=3300000 #давление открытия аварийного клапана, Па

#геометрические параметры:
V_ballon=10. #объем баллона, м3
d_valve1=0.01 #диаметр отверстия 1го клапана, м
d_valve2=0.009 #диаметр отверстия 2го клапана, м
d_valve3=0.011 #диаметр отверстия 3го клапана, м
d_safety_valve=0.03 #диаметр аварийного клапана, м

S1=np.pi*(d_valve1/2)**2#площадь сечения клапана
S2=np.pi*(d_valve2/2)**2
S3=np.pi*(d_valve3/2)**2
Ssafety=np.pi*(d_safety_valve/2)**2

#временнЫе настройки расчета
T=100 #общее время расчета
dt=0.001 #шаг расчета по времени
time_compressor= {0: True,#время включения/выключения компрессора
                 1:True,
                 45:False,
                 55:True,
                90:False}

time_valve1= {0: False,#время открытия/закрытия клапана (False - закрыт, True- открыт)
                 30:True,
                85:False}

time_valve2= {0: False,#время открытия/закрытия клапана (False - закрыт, True- открыт)
                 60:True,
                90:False}

time_valve3= {0: False,#время открытия/закрытия клапана (False - закрыт, True- открыт)
                 80:True,
                95:False}

#первоначальные приближения
P_ballon0=101325 #давление в баллоне
T_ballon0=288.15
compr_status=True
v1_status=False
v2_status=False
v3_status=False
G_safety_valve_status=False
G_v1=0#расходы через клапаны
G_v2=0
G_v3=0
G_safety_valve=0
P_ballon_current=P_ballon0
T_ballon_current=T_ballon0
Rho_ballon_current=P_ballon_current/R/T_ballon_current
M_current=Rho_ballon_current*V_ballon #масса воздуха в баллоне
lam_previous=0.0001

#Массивы для записи результатов:
time=[0]
res_P_ballon=[P_ballon0]
res_T_ballon=[T_ballon0]
res_Rho_ballon=[Rho_ballon_current]
res_G_in=[0]
res_G_v1=[G_v1]
res_G_v2=[G_v2]
res_G_v3=[G_v3]
res_G_safety_valve=[G_safety_valve]
res_safety_valve=[False]
res_M=[np.nan]

#газодинамическая функция отношение статического к плоному давлению в зависимости от безразмерной скорости
def Pi(lamda,k):
    return (1-lamda**2*((k-1)/(k+1)))**(k/(k-1))

def make_plots():
    fig3=mc.Chart(points_for_plot=[{'x':time,'y':res_G_in,'label':'G_compr','ls':'dashed','c':'gray'},
                                   {'x':time,'y':res_G_v1,'label':'G_v1','c':'green'},
                                   {'x':time,'y':res_G_v2,'label':'G_v2','c':'blue'},
                                   {'x':time,'y':res_G_v3,'label':'G_v3','c':'red'},
                                   {'x':time,'y':res_G_safety_valve,'label':'G_safety_valve','lw':0.1,'ls':'dashed','c':'black'}],xlabel='time, s',ylabel='G, kg/s', title='Massflows', dpi=150,figure_size=(5,5))
    fig4=mc.Chart(points_for_plot=[{'x':time,'y':res_P_ballon,'label':'P_ballon','c':'red'},],xlabel='time, s',ylabel='P, Pa', title='Pressure', dpi=150,figure_size=(5,5))
    fig5=mc.Chart(points_for_plot=[{'x':time,'y':res_T_ballon,'label':'T_ballon','c':'red'},],xlabel='time, s',ylabel='T, K', title='Temperature', dpi=150,figure_size=(5,5))
    fig6=mc.Chart(points_for_plot=[{'x':time,'y':res_Rho_ballon,'label':'Rho_ballon','c':'red'},],xlabel='time, s',ylabel='Rho, kg/m3', title='Density', dpi=150,figure_size=(5,5))
    fig7=mc.Chart(points_for_plot=[{'x':time,'y':res_M,'label':'M_ballon','c':'red'},],xlabel='time, s',ylabel='Mass of air in ballon, kg', title='mass', dpi=150,figure_size=(5,5))


#Расчет
try:
    for i in range(1,int(T/dt)):
        t=i*dt
        time.append(t)
        if i%1000==0:
            print(f"step calculated: {i}")

        # проверяем статус работы компрессора
        for t_, compr_status_ in time_compressor.items():
            if (abs(t_ - t) < dt * 0.5):
                compr_status = compr_status_

        # проверяем статус работы клапанов 1-3
        for t_, v1_status_ in time_valve1.items():
            if (abs(t_ - t) < dt * 0.5):
                v1_status = v1_status_
        for t_, v2_status_ in time_valve2.items():
            if (abs(t_ - t) < dt * 0.5):
                v2_status = v2_status_
        for t_, v3_status_ in time_valve3.items():
            if (abs(t_ - t) < dt * 0.5):
                v3_status = v3_status_

        #расход воздуха из компрессора в баллон
        if compr_status:
            PR= P_ballon_current / P_atm
            G_in=float(compressor_curve_for_calculation(PR))
        else:
            G_in=0.

        #считаем параметры газа на выходе из компрессора, т.е. то, что попадает в баллон
        T_out_is=T_compr_in*(P_ballon_current / P_atm)**((k-1)/k)#из уравнения адиабатического процесса идеальная температура без учета кпд
        L_is=Cp*(T_out_is-T_compr_in)#изоэнтропическая работа компрессора
        L_real=L_is/eff_compr #реальная работа компрессора
        T_out=L_real/Cp+T_compr_in # реальная темпеартура на выходе из компрессора

        #найдем расходы через клапаны 1-3:
        #сначала из газодинамических функций найдем безразмерную скорость потока
        if P_atm/P_ballon_current<0.6: #т.к. заведомо известно, что для воздуха критический перепад давления около 0,55-0,6, то логически отсеим вариант с бОльшим перепадом
            lam=1.
        else:
            f = lambda lam: Pi(lam,k) - P_atm / P_ballon_current #вспомогательная функция для поиска скорости истечения на основе известного перепада давления
            # res = root_scalar(f, x0=lam_previous, x1=lam_previous*1.0001, method='secant') #обычно чуть более устойчивый но медленный метод
            res = root_scalar(f, bracket=[0.,2.4], method='toms748') #считает быстрее, но нужно задавать границы поиска корня - можно придумать как их задавать, чтобы обеспечить гарантиронное решение, в данном случае я просто задал границы от 0 до 2
            lam=res.root #безразмерная скорость при истечении из клапанов
            if lam>1:
                lam=1

        lam_previous=lam

        #далее для каждого клапана ищем скорости
        a_critical=(2*k/(k+1)*R*T_ballon_current)**0.5
        vel=lam*a_critical
        if v1_status:
            G_v1 = Rho_ballon_current * vel * S1
        else:
            G_v1=0.

        if v2_status:
            G_v2 = Rho_ballon_current * vel * S2
        else:
            G_v2=0.

        if v3_status:
            G_v3 = Rho_ballon_current * vel * S3
        else:
            G_v2=0.

        #проверка открыт ли аварийный клапан:
        if P_ballon_current>P_safety_valve:
            G_safety_valve=Rho_ballon_current * vel * Ssafety
            G_safety_valve_status=True
        elif P_ballon_current<P_safety_valve*0.99: #сделал, что-то вроде гистерезиса, чтобы клапан обратно закрывался не сразу, а с запаздыванием - так графики выглядят более "по-человечески"
            G_safety_valve=0.
            G_safety_valve_status = False

        #далее считаем параметры в баллоне после смешения новой порции воздуха из компрессора с воздухом в баллоне:
        M_compr_in=G_in*dt
        M_out=(G_v1+G_v2+G_v3+G_safety_valve)*dt
        dM=M_compr_in-M_out
        M_new=M_current+dM #новое значение массы с учетом подводов и отборов

        H_ballon=Cp*T_ballon_current*M_current #энтальпия в баллоне в начале текущего шага
        H_out=Cp*T_ballon_current*M_out #энтальпия выводимая из баллона с воздухом через клапана
        H_in=Cp*T_out*M_compr_in #энтальпия подводимая в баллон от компрессора
        H_new=H_ballon+H_in-H_out #энтальпия в баллоне в конце текущего шага

        #далее условно считаем адиабатический процесс изменения параметров в баллоне с учетом новой массы
        # T_new=(T_out*M_compr_in+T_ballon_current*M_new)/M_new
        T_new=H_new/Cp/M_new
        Rho_new=M_new/V_ballon
        P_new=R*T_new*Rho_new
        # P_new=P_ballon_current*(T_new/T_ballon_current)**(k/(k-1))
        # Rho_new=P_new/R/T_new

        P_ballon_current=P_new
        T_ballon_current=T_new
        Rho_ballon_current=Rho_new
        M_current=M_new

        #сохраняем результаты:
        res_P_ballon.append(P_new)
        res_T_ballon.append(T_new)
        res_Rho_ballon.append(Rho_new)
        res_G_in.append(G_in)
        res_G_v1.append(G_v1)
        res_G_v2.append(G_v2)
        res_G_v3.append(G_v3)
        res_G_safety_valve.append(G_safety_valve)
        res_safety_valve.append(G_safety_valve_status)
        res_M.append(M_new)


except:
    pass
    # make_plots()
    # plt.show()

make_plots()


#на всяеий случай проверим газодинамическую функцию по давлению:
# x_array=np.linspace(0,3.,100)
# y_array=[Pi(lam,k) for lam in x_array]
# fig2=mc.Chart(points_for_plot=[{'x':x_array,'y':y_array,'label':'Pi','ls':'dashed','c':'red'}],xlabel='lamda',ylabel='Pi', dpi=150,figure_size=(5,5))

plt.show()