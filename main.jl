using Plots
gr()
using LinearAlgebra

S0 = 1/2
L0 = 3
j0 = abs(S0-L0)

println("j0 = $j0")

gj = 1 + (j0*(j0+1) + S0*(S0+1) - L0*(L0+1))/(2*j0*(j0+1))

j0 = 2*j0 + 1

println("j0 = $j0")

E = zeros(Int(j0),Int(j0))
JZ = zeros(Int(j0),Int(j0))
JP = zeros(Int(j0),Int(j0))
JM = zeros(Int(j0),Int(j0))

answ = zeros(100,100)
answ = reshape(Complex.(answ),100,100)

display(answ)

for x in 1:100
    T(x) = 0.01 + 100*(x-100)/100
    Ha = [0 0 1]
    Hmag = .1 #Matlabの .1 は Juliaの 0.1と同じ？
    mem = T(x)
    Hmag = 10000 * gj * 0.92740154*10^(-20) / 1.380658 * 10^(-16) * Hmag
    
    #println("clear01")

    for k in 1:Int(j0)
        for j in 1:Int(j0)
            jnum = (j-1)/2
            mi = k - jnum -1
            mj = j - jnum -1
            if k == j
                E[k,j] = 1
                JZ[k,j] = mi
            end
            if k == j-1
                JP[k,j] = ((jnum + mi + 1) * (jnum - mi))^(0.5)
            end
            if k == j + 1
                JM[k,j] = ((jnum + mi) * (jnum - mi + 1))^(0.5)
            end
        end
    end
    
    #println("clear02")

    #@time display(E)
    #@time display(JZ)
    #@time display(JP)
    #@time display(JM)

    JX = 0.5 * (JP + JM)
    JY = (-0.5*im)*(JP - JM)
    O40 = 35 * JZ^4 + (-30 * jnum * (jnum + 1) + 25) * JZ^2 + (-6 * jnum * (jnum + 1) + 3 * jnum^2 * (jnum + 1)^2) * E
    O44 = 0.5 * (JP^4 + JM^4)
    O4 = O40 + 5 * O44
    HZeeman = -Hmag * JX
    
    #println("clear03")

    B4 = 1 #Γ7 基底状態
    #B4 = -1 #Γ8 基底状態

    Hcrys = 60 * B4 * (O40 + 5 * O44)
    Hami = HZeeman + Hcrys
    
    #println("clear04")

    #固有値 XR と ユニタリ行列D
    XR,D = eigen(Hami)

    XI = XR'

    dis = zeros(6,6)
    
    #println("clear05")

    for k in 1:6
        for j in 1:6
            if k == j
                dis[k,j] = exp(-D[k,k]/T(x))
            else
                dis[k,j] = 0
            end
        end
    end
    
    #println("clear06")

    EIG1 = D[1,1]
    EIG2 = D[2,2]
    Eav = tr(D*dis)/tr(dis)
    
    #println("clear07")

    Dum = XI * JX * XR
    JXav = tr(Dum * dis)/tr(dis)
    
    #println("clear08")

    Dum = XI * JY * XR
    JYav = tr(Dum * dis)/tr(dis)
    
    #println("clear09")

    Dum = XI * JZ * XR
    JZav = tr(Dum * dis)/tr(dis)
    
    #println("clear09")

    J001 = JXav
    J110 = (JXav + JYav) / sqrt(2)
    
    #println("clear010")

    chi = (JXav + JYav)/(Hmag + eps(1.0)) #epsの使い方は疑問
    
    #println("clear011")

    invchi = 1/chi
    
    #println("invchi = $invchi")
    
    #println("clear012")

    answ[x,1] = jnum #jnum ?
    answ[x,2] = JXav
    answ[x,3] = EIG1
    answ[x,4] = EIG2
    answ[x,5] = Eav
    answ[x,6] = J001
    answ[x,7] = invchi
    
    #println("clear013")
end

display(answ)

display(real(answ[:,1]))
display(real(answ[:,2]))
display(real(answ[:,3]))


plot(real(answ[:,1]),real(answ[:,2]))
xlabel!("H(T)")
ylabel!("M")

scatter(real(answ[:,1]),real(answ[:,3]))
#plot!(real(answ[:,1]),real(answ[:,4]))
xlabel!("H(T)")
ylabel!("M")