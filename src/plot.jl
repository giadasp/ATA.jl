function pgfplot()
 Plots.pgfplotsx()
end

function PlotATA!(ATAmodel::model, IIFf, ICFf, design::Matrix{Float64}; simPool=Float64[], results_folder="RESULTS")
    T=ATAmodel.Settings.T
    ThetasPlot=collect(range(-4,stop=4,length=101))
    if ATAmodel.Settings.nGroups>1
        n_colors=Plots.Colors.range(Plots.Colors.colorant"red",Plots.Colors.colorant"green", length = ATAmodel.Settings.nGroups)
    else
        n_colors = [Plots.Colors.colorant"blue"]
    end
    colors=vcat([[n_colors[g] for t=1:ATAmodel.Settings.Tg[g]] for g=1:ATAmodel.Settings.nGroups]...)
    Plots.plot(ThetasPlot,IIFf',
    xlims=(-4,4),
    xticks = -4:1:4,
    titlefontsize=16,
    size=(500,400),
    yticks=0:2:maximum(IIFf)+2,
    ylims=(0,maximum(IIFf)+2),
    tickfontsize=14,
    guidefontsize=14,
    legendfontsize=16,
    xtickfontrotation=45,
    thickness_scaling=0.8,
    foreground_color_border=:black,
    palette=colors,
    linewidth=1.0,
    alpha=0.5,
    labels=permutedims([string("t",t) for t=1:T]),
    );
    Plots.yaxis!(L"TIF({\theta})");
    Plots.xaxis!(L"{\theta}");
    Plots.savefig(string(results_folder,"/TIFPlot.pdf"))
    maxscore=[sum(design[:,t]) for t=1:T]
    ICFf=[ICFf[t,i]./maxscore[t] for t=1:T, i=1:101]
    Plots.plot(ThetasPlot,ICFf',
    xlims=(-4,4),
    xticks = -4:1:4,
    size=(500,400),
    yticks=0:0.1:1,
    ylims=(0,1),
    titlefontsize=16,
    tickfontsize=14,
    guidefontsize=14,
    legendfontsize=16,
    xtickfontrotation=45,
    thickness_scaling=0.8,
    foreground_color_border=:black,
    palette=colors,
    linewidth=1.0,
    alpha=0.5,
    labels=permutedims([string("t",t) for t=1:T]),
    );
    Plots.yaxis!(L"ICF({\theta})");
    Plots.xaxis!(L"{\theta}");
    Plots.savefig(string(results_folder,"/ICFPlot.pdf"))
    if size(simPool,1)>0
        IIFtrue=Vector{Vector{Float64}}(undef,T)
        for t=1:T
            IIFtrue[t]=ItemInfoFun(simPool,ThetasPlot, model = ATAmodel.Settings.IRT.model, parametrization = ATAmodel.Settings.IRT.parametrization, D = ATAmodel.Settings.IRT.D)'*design[:,t]
        end
    end
end
function PlotATA_CC!(ATAmodel::model, IIFf, ICFf, design::Matrix{Float64}; simPool=Float64[], results_folder="RESULTS")
    T=ATAmodel.Settings.T
    alphaR=Int(ceil(ATAmodel.Obj.AuxInt*(ATAmodel.Obj.AuxFloat)))
    IIF_plot=Vector{Array{Float64,3}}(undef,T)
    ICF_CC_plot=Vector{Array{Float64,3}}(undef,T)
    ThetasPlot=collect(range(-4,stop=4,length=101)) #nqp values in interval/r/n",
    if isfile("BSPar.jld2")
        JLD2.@load "BSPar.jld2" BSPar
    else
        error(string("file BSPar.jld2 not found"))
    end
    BSa=Matrix(BSPar[2])[:,2:end]
    BSb=Matrix(BSPar[1])[:,2:end]
    R=ATAmodel.Obj.AuxInt
    for t=1:T
        println(t)
        IIF_plot[t]=zeros(101,ATAmodel.Settings.nItems,R)
        ICF_CC_plot[t]=zeros(101,ATAmodel.Settings.nItems,R)
        for r=1:R
            if ATAmodel.Settings.IRT.model=="1PL"
                df=DataFrame(b=BSb[:,r]) #nqp values in interval\r\n",
            elseif ATAmodel.Settings.IRT.model=="2PL"
                df=DataFrame(a=BSa[:,r],b=BSb[:,r]) #nqp values in interval\r\n",
            elseif ATAmodel.Settings.IRT.model=="3PL"
                df=DataFrame(a=BSa[:,r],b=BSb[:,r],c=BSc[:,r])
            end
            for k=1:101
                IIF_plot[t][k,:,r]=ItemInfoFun(df,ThetasPlot[k],model=ATAmodel.Settings.IRT.model, parametrization=ATAmodel.Settings.IRT.parametrization, D=ATAmodel.Settings.IRT.D) # IxK[t]
                ICF_CC_plot[t][k,:,r]=ItemCharFun(df,ThetasPlot[k],model=ATAmodel.Settings.IRT.model, parametrization=ATAmodel.Settings.IRT.parametrization, D=ATAmodel.Settings.IRT.D)# IxK[t]
            end
        end
    end
    #TIF=Array{Array{Float64,2},1}(undef,T)
    IIFdesigntoplot=Array{Array{Float64,2},1}(undef,T)
    for t in 1:T
        TIF=Array{Float64,2}(undef,101,R)
        IIFdesigntoplot[t]=Array{Float64,2}(undef,101,6)
        for k=1:101
            for r in 1:R
                TIF[k,r]=(IIF_plot[t][k,:,r]'*design[:,t])#,[0,0.25,0.5,0.75,1,α])[1:6]
            end
            IIFdesigntoplot[t][k,:]=sort(TIF[k,:])[[1,Int(ceil(ATAmodel.Obj.AuxInt*0.25)),Int(ceil(ATAmodel.Obj.AuxInt*0.5)),Int(ceil(ATAmodel.Obj.AuxInt*0.75)),ATAmodel.Obj.AuxInt,alphaR]]
        end
        DelimitedFiles.writedlm(string(results_folder,"/IIFdesigntoplot_",t,".csv"),IIFdesigntoplot[t])
    end
    for t in 1:T
        #plot(IIFdesigntoplot[t][:,1],IIFdesigntoplot[t][:,2],seriestype=:scatter)
        Plots.plot(size=(500,400),yticks=0.0:2.0:maximum(IIFf)+2,ylims=(0,maximum(IIFf)+2),xlims=(-4.0,4.0));
        if size(simPool,1)>0
            Plots.plot!(ThetasPlot,IIFtrue[t],size=(500,400),yticks=0:2:maximum(IIFf)+2,ylims=(0,maximum(IIFf)+2),tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,linestyle = :solid,linecolor=:violetred4,label="True");
        end
        Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,5],
        titlefontsize=16, size=(500,400),tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,
        linestyle = :dot,linecolor=:darkcyan,label="Max");
        Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,4],titlefontsize=16,tickfontsize=14,guidefontsize=14,legendfontsize=16, xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,yticks=0:2:maximum(IIFf)+2,ylims=(0,maximum(IIFf)+2),linestyle = :dash,linecolor=:darkcyan,label="75-Qle");
        Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,3],titlefontsize=16,tickfontsize=14,guidefontsize=14,legendfontsize=16,  xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,linestyle = :solid,linecolor=:darkcyan,label="Median");
        Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,2],titlefontsize=16,tickfontsize=14,guidefontsize=14,legendfontsize=16,  xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,linestyle = :dash,linecolor=:darkcyan,label="25-Qle");
        Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,6],titlefontsize=16,tickfontsize=14,guidefontsize=14,legendfontsize=16,  xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,linestyle = :dashdotdot,linecolor=:indigo,label=L"{\alpha}-Qle");
        Plots.plot!(ThetasPlot,IIFdesigntoplot[t][:,1],titlefontsize=16,tickfontsize=14,guidefontsize=14,legendfontsize=16,  xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,linestyle = :dot,linecolor=:darkcyan,label="Min");
        Plots.yaxis!(L"TIF({\theta})");
        Plots.xaxis!(L"{\theta}");
        Plots.plot!(ThetasPlot,IIFf[t,:],tickfontsize=12, markersize=4, xtickfontrotation=45, thickness_scaling=0.8,foreground_color_border=:black, linewidth=1.0,colour=[:black],label="estimated");
        Plots.savefig(string(results_folder,"/",t,"_TIFPlot.pdf"))
    end
end
