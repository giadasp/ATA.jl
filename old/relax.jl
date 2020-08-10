function minF(ATAmodel::ATA.model,IIF,OptFeas::Float64,opt::NLopt.Opt; findFeasible=true)
    FScounts=ATAmodel.Settings.FS.Counts*ones(Float64,ATAmodel.Settings.T)'
    function myf(x::Vector,grad::Vector)
        infeas=zeros(Float64,ATAmodel.Settings.T)
        obj=zeros(Float64,ATAmodel.Settings.T)
        ol=zeros(Float64,ATAmodel.Settings.T)
        x=reshape(x,ATAmodel.Settings.nFS,ATAmodel.Settings.T)
        x=round.(x)
        for v=1:ATAmodel.Settings.T
            infeas[v], x_I_v=ATA.checkFeas(ATAmodel.Settings.FS,ATAmodel.Constraints[v],x[:,v],ATAmodel.Settings.nFS,ATAmodel.Settings.nItems,v)
            # if findFeasible==false
            #     if ATAmodel.Settings.OptType=="MAXIMIN"
            #         obj[v]=ATA.evalTIFMMv(x_I_v,IIF[v])
            #     elseif ATAmodel.Settings.OptType=="CC"
            #         obj[v]=ATA.evalTIFCCv(x_I_v,IIF[v];α=ATAmodel.Obj.AuxFloat)
            #     end
            # end
        end
        ol=ATA.evalOverlap(x,FScounts,ATAmodel.Settings.olMax,ATAmodel.Settings.T,ol)
        iu=sum(x,dims=2)-ATAmodel.IU.Max
        iu=iu[iu.>0]
        if size(iu,1)==0
            iu=0
        else
            iu=sum(iu)
        end
        if length(grad)>0
            for i=1:ATAmodel.Settings.nFS*ATAmodel.Settings.T
                grad[i]=0.0
            end
        end
        return (1-OptFeas)*(sum(infeas+ol)+iu)-OptFeas*minimum(obj)
    end
    opt.min_objective = myf
    opt_f = Array{Cdouble}(undef,1)
    x=ones(Float64,ATAmodel.Settings.nFS*ATAmodel.Settings.T).*(1/2)
    ccall((:nlopt_optimize,NLopt.libnlopt), NLopt.Result, (NLopt._Opt, Ptr{Cdouble}, Ptr{Cdouble}), opt, x, opt_f)
    return x
end

function relaxATA(ATAmodel::ATA.model,IIF,OptFeas::Float64)
    opt=NLopt.Opt(:LN_COBYLA,ATAmodel.Settings.nFS*ATAmodel.Settings.T)
    opt.lower_bounds = zeros(ATAmodel.Settings.nFS*ATAmodel.Settings.T)
    opt.upper_bounds = ones(ATAmodel.Settings.nFS*ATAmodel.Settings.T)
    opt.xtol_rel = 1e-5
    opt.maxtime = 10.0
    opt.ftol_rel=  1e-6
    x=maxF(ATAmodel,IIF,OptFeas,opt)
    return x
end
