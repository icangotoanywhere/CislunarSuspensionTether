import copy
import numpy as np
import math
import wlph

N=9

x_ms = 3.7826E8 #月面距离地月公共质心的距离（月心减去月球半径），单位m
x_end= 0.4E8 #天梯端点所处的位置，单位m
l_s=(x_ms-x_end)/N #每段天梯的长度，单位m


sizepop=600    #种群规模
maxiter=50000    #最大繁殖代数

tSize=4        #锦标赛规模

lenchrom=N+1    #质点个数 

#以下子函数调用etabs进行频率分析

def Analysis1(X,Y):

    fobj=np.sum(X)+np.sum(Y)


    return fobj




#以下子程序，实现种群变异（基因层面）
def mutation(chrom,sizepop,lenchrom):
    chrom1=[]
    ffobj=[]
    for i in range(sizepop):
        X=chrom[i][0][:]
        Y=chrom[i][1][:]
        rr=np.random.rand(1)
        k1=math.floor((lenchrom-1)*rr[0]) #确定变异位置
        for k in range(k1,lenchrom):        
            ddd=0
            while ddd<l_s**2 or ddd>1.21*l_s**2:
                rad=np.random.rand(2)

                XN=X[k-1]+(2.0*rad[0]-1)*l_s
                YN=Y[k-1]+(2.0*rad[1]-1)*l_s
                ddd=(XN-X[k-1])**2+(YN-Y[k-1])**2
            X[k]=XN
            Y[k]=YN


        x0=X[1:]
        y0=Y[1:]
        fcrm=wlph.F(x0,y0)


        chrom1.append([X,Y])


        ffobj.append(fcrm)
    chrom=chrom1

    return chrom,ffobj



#以下语句，锦标赛选择
def tournament(chrom,ffobj,sizepop,tSize):
    chrom2=[]
    chrom2.append([Xopt,Yopt])

    for i in range(1,sizepop):
        ff0=1E15
     
        for j in range(tSize):
            rr=np.random.rand(1)[0]
            idr=math.floor(rr*sizepop)

            
            if ffobj[idr]<ff0:
                ff0=ffobj[idr]
                idrc=idr

        chrom2.append(chrom[idrc])
    chrom=chrom2
    
    return chrom


#以下语句，生成初始种群,计算初始种群适应度

chrom=[]
ffobj=[]

for j in range(sizepop):
    X=[0]
    Y=[0]
    for i in range(1,lenchrom):

        ddd=0
        while ddd<l_s**2 or ddd>1.21*l_s**2:
            rad=np.random.rand(2)

            XN=X[i-1]+(2.0*rad[0]-1)*l_s
            YN=Y[i-1]+(2.0*rad[1]-1)*l_s
            ddd=(XN-X[i-1])**2+(YN-Y[i-1])**2


        X.append(XN)
        Y.append(YN)

    x0=X[1:]
    y0=Y[1:]
    fcrm=wlph.F(x0,y0)

    chrom.append([X,Y])
    ffobj.append(fcrm)

print('ffobj',ffobj)

   

#以下语句，选出适应度值和最优个体

ffobjopt=min(ffobj)

print('初始种群ffobjopt',ffobjopt)

Xopt=[]
Yopt=[]
for i in range(lenchrom):
    Xopt.append(chrom[ffobj.index(min(ffobj))][0][i])
    Yopt.append(chrom[ffobj.index(min(ffobj))][1][i])


print('初始种群最优解',Xopt,Yopt)





#以下为主循环程序
for i in range(maxiter):
    print("第",i,"次循环开始")

    chrom=tournament(chrom,ffobj,sizepop,tSize)

    chrom,ffobj=mutation(chrom,sizepop,lenchrom)


    if ffobjopt>min(ffobj):
        ffobjopt=min(ffobj)
        Xopt=chrom[ffobj.index(min(ffobj))][0]
        Yopt=chrom[ffobj.index(min(ffobj))][1]
        print("最优解是")
        print(Xopt,Yopt)

    if ffobjopt<min(ffobj):
        chrom[0]=[Xopt,Yopt]
        ffobj[0]=ffobjopt

          
    print("以上是第",i,"次循环")
    print("最优解是")
    print(Xopt,Yopt)
    print("ffobjopt=",ffobjopt)


ffinal=Analysis1(Xopt,Yopt)

print("全部结束！累死哥了，差点跑掉一个电容！")

print("ffobjopt=",ffobjopt)
