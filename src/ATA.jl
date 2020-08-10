__precompile__()
module ATA

import Printf
import CSV
import DelimitedFiles
import DataFrames
import Distributions
import Distributed
import Random
import LinearAlgebra
import JLD2
import StatsBase
import LaTeXStrings.@L_str
import Base.Threads
import Dates
import MathOptInterface
import JuMP
import Dash
import DashHtmlComponents
import DashCoreComponents
import Plots
import PGFPlotsX
using Requires

const DataFrame = DataFrames.DataFrame
const Distribution = Distributions.Distribution
const readdlm = DelimitedFiles.readdlm
const writedlm = DelimitedFiles.writedlm
const dot = LinearAlgebra.dot
const sample = StatsBase.sample
const quantile = StatsBase.quantile

function mean(x::Vector{Real})
	n = size(x, 1)
	if n > 0
		return sum(x) / n
	else
		return 0
	end
end

mutable struct InputSettings
	T::Vector{Int64}
	nItems::Int64
	nGroups::Int64
	Groups::Vector{String}
	IRTmodel::String
	IRTparameters::Vector{String}
	IRTparametrization::String
	IRTD::Float64
	EnemySetsVar::Vector{String}
	FriendSetsVar::Vector{String}
	ItemUse_min::Vector{Int64}
	ItemUse_max::Vector{Int64}
	length_min::Vector{Int64}
	length_max::Vector{Int64}
	length_weight::Vector{Float64}
	ExSVar::Vector{String}
	ExSPts::Vector{Vector{Float64}}
	ExS_min::Vector{Vector{Float64}}
	ExS_max::Vector{Vector{Float64}}
	meanVars::Vector{Vector{String}}
	meanVars_min::Vector{Vector{Float64}}
	meanVars_max::Vector{Vector{Float64}}
	sumVars::Vector{Vector{String}}
	sumVars_min::Vector{Vector{Float64}}
	sumVars_max::Vector{Vector{Float64}}
	OptType::String
	OptPts::Vector{Vector{Float64}}
	AuxInt::Int64
	AuxFloat::Float64
	CATEGORIES::Vector{String}
	InputSettings(T, nItems, nGroups, Groups,
	IRTmodel, IRTparameters, IRTparametrization, IRTD,
	EnemySetsVar, FriendSetsVar,
	ItemUse_min, ItemUse_max,
	length_min, length_max, length_weight,
	ExSVar, ExSPts, ExS_min,ExS_max,
	meanVars, meanVars_min, meanVars_max,
	sumVars, sumVars_min, sumVars_max,
	OptType, OptPts, AuxInt, AuxFloat,
	CATEGORIES) = new(T, nItems, nGroups, Groups,
	IRTmodel, IRTparameters, IRTparametrization, IRTD,
	EnemySetsVar, FriendSetsVar,
	ItemUse_min, ItemUse_max,
	length_min, length_max, length_weight,
	ExSVar, ExSPts, ExS_min, ExS_max,
	meanVars, meanVars_min, meanVars_max,
	sumVars, sumVars_min, sumVars_max,
	OptType, OptPts, AuxInt, AuxFloat,
	CATEGORIES)
	InputSettings() = new(Int64[], zero(Int64), zero(Int64), String[],
	"", String[], "at-b", 1.0,
	String[], String[],
	Int64[], Int64[],
	Int64[], Int64[], Float64[],
	String[], Float64[], Float64[], Float64[],
	Vector{Vector{String}}(undef, 0), Vector{Vector{Float64}}(undef, 0), Vector{Vector{Float64}}(undef, 0),
	Vector{Vector{String}}(undef, 0), Vector{Vector{Float64}}(undef, 0), Vector{Vector{Float64}}(undef, 0),
	"MAXIMIN",Vector{Vector{Float64}}(undef, 0), zero(Int64), zero(Float64),
	String[])
end

mutable struct IRT
	model::String
	parameters::DataFrame
	parametrization::String
	D::Float64
	metric::Vector{Float64} # [meanTarget,stdTarget]
	X::Vector{Float64}
	W::Vector{Float64}
	IRT() = new("2PL", DataFrame(), "at-b", 1.0, [0.0,1.0], zeros(1), zeros(1))
	IRT(model, parameters, parametrization, D, metric, X, W) = new(model, parameters, parametrization, D, metric, X, W)
end

mutable struct FS
	Var::Vector{Symbol}
	Sets::Vector{String}
	Counts::Vector{Int64}
	Items::Vector{Vector{Int64}}
	FS(Var,Sets,Counts,Items) = new(Var,Sets,Counts,Items)
	FS()=new(Symbol[],String[],Int64[],Vector{Vector{Int64}}(undef,0))
end

mutable struct ES
	Var::Vector{Symbol}
	Names::Vector{String}
	Sets::Vector{Vector{Int64}}
	ES(Var,Names,Sets) = new(Var,Names,Sets)
	ES()=new(Symbol[],String[],Vector{Vector{Int64}}(undef,0))
end

mutable struct ExS
	Var::Symbol
	Val::Matrix{Float64}
	Min::Vector{Float64}
	Max::Vector{Float64}
	Pts::Vector{Float64}
	ExS(Var, Val, Min, Max, Pts) = new(Var, Val, Min, Max, Pts)
	ExS() = new(Symbol(""), zeros(Float64,0 , 0), Float64[], Float64[], Float64[])
end

mutable struct Settings
	nItems::Int64
	nFS::Int64
	Bank::DataFrame
	IRT::IRT
	ThetaBounds::Vector{Vector{Float64}}
	forced0::Vector{Vector{Bool}}
	nGroups::Int64
	T::Int64
	Tg::Vector{Int64}
	OptType::String
	FS::FS #friend Settings
	ES::ES #enemy Settings
	olMax::Matrix{Float64}
	Settings(nItems, nFS, Bank, IRT, ThetaBounds, forced0, nGroups, T, Tg, OptType, FS, ES, olMax) = new(nItems, nFS, Bank, IRT, ThetaBounds, forced0, nGroups, T, Tg, OptType, FS, ES, olMax) #no pattern mode
	Settings() = new(0, 0, DataFrame(), IRT(), [[-6.0,6.0]], Vector{Vector{Bool}}(undef, 0), 1, 1, [1], "MAXIMIN", FS(), ES(), zeros(Int64, 0, 0))
end

mutable struct Neighbourhood
	x::Matrix{Float64}
	f::Float64
	obj::Vector{Float64}
	infeas::Vector{Float64}
	ol::Vector{Float64}
	iu::Float64
	Neighbourhood(x, f, obj, infeas, ol, iu)=new(x, f, obj, infeas, ol, iu)
	Neighbourhood() = new(Matrix{Float64}(undef, 0, 0), Inf, Float64[], Float64[], Float64[], zero(Float64))
end

mutable struct Constraint
	length_min::Int64
	length_max::Int64
	ExS::ExS
	meanVars::Vector{Symbol} #constrain the mean to be
	meanVars_min::Vector{Float64}
	meanVars_max::Vector{Float64}
	sumVars::Vector{Symbol}#constrain the sum to be
	sumVars_min::Vector{Float64}
	sumVars_max::Vector{Float64}
	catConstrA::Matrix{Float64}
	catConstrb::Vector{Float64}
	#olMax::Matrix{Int64}
	Constraint(length_min, length_max, ExS, meanVars, meanVars_min, meanVars_max, sumVars, sumVars_min, sumVars_max, catConstrA, catConstrb) = new(length_min, length_max, ExS, meanVars, meanVars_min, meanVars_max, sumVars, sumVars_min, sumVars_max, catConstrA, catConstrb)
	Constraint() = new(zero(Int64), 10000000, ExS(), Vector{Symbol}(undef, 0), Vector{Float64}(undef, 0), Vector{Float64}(undef, 0), Vector{Symbol}(undef, 0), Vector{Float64}(undef, 0), Vector{Float64}(undef, 0), Matrix{Float64}(undef, 0, 0), Vector{Float64}(undef, 0))
end

mutable struct Obj
	Sense::String
	OptPts::Vector{Vector{Float64}}
	AuxInt::Int64
	AuxFloat::Float64
	Obj(Sense, OptPts, AuxInt, AuxFloat) = new(Sense, OptPts, AuxInt, AuxFloat)
	Obj() = new("Max", Vector{Vector{Float64}}(undef, 0), zero(Int64), zero(Float64))
end

mutable struct Output
	Categories::Vector{Symbol}
	Quantitative::Vector{Symbol}
	SummQuan::Vector{Function}
	Design::Matrix{Float64}
	f::Float64
	Feas::Vector{Float64}
	TIF::Vector{Float64}
	ElapsedTime::Float64
	Neighbourhoods::Vector{Neighbourhood}
	Output(Categories, Quantitative, summQuan, Design, f, Feas, TIF, ElapsedTime, Neighbourhoods)=new(Categories, Quantitative, SummQuan, Design, f, Feas, TIF, ElapsedTime, Neighbourhoods)
	Output() = new(Vector{Vector{Symbol}}(undef, 0), Vector{Vector{Symbol}}(undef, 0), Function[], zeros(Float64, 0, 0), zero(Float64), Float64[], Float64[], zero(Float64), Neighbourhood[])
end

mutable struct IU
	Min::Vector{Int64}
	Max::Vector{Int64}
	IU(Min, Max) = new(Min, Max)
	IU() = new(Int64[], Int64[])
end

mutable struct model
	Settings::Settings
	Constraints::Vector{Constraint}
	IU::IU
	Obj::Obj
	Output::Output
	model(Settings, Constraints, IU, Obj, Output)=new(Settings, Constraints, IU, Obj, Output)
	model() = new(Settings(), Constraint[], IU(), Obj(), Output())
end

mutable struct opt
	x::Matrix{Float64}
	f::Float64
	obj::Vector{Float64}
	infeas::Vector{Float64}
	ol::Vector{Float64}
	iu::Vector{Float64}
	opt(x, f, obj, infeas, ol, iu)=new(x, f, obj, infeas, ol, iu)
	opt() = new(Matrix{Float64}(undef, 0, 0), Inf, Float64[], Float64[], Float64[], Float64[])
end

include("utils.jl")
include("build.jl")
include("opt.jl")
include("plot.jl")
include("out.jl")
include("app.jl")
function __init__()
	@require Cbc = "9961bab8-2fa3-5c5a-9d89-47fab24efd76" begin
		function add_Cbc!(m)
			JuMP.set_optimizer(m, Cbc.Optimizer)
		end
	end
	@require Gurobi = "2e9cd046-0924-5485-92f1-d5272153d98b" begin
		function add_Gurobi!(m)
			JuMP.set_optimizer(m, Gurobi.Optimizer)
		end
	end
	@require CPLEX = "a076750e-1247-5638-91d2-ce28b192dca0" begin
		function add_CPLEX!(m)
			JuMP.set_optimizer(m, CPLEX.Optimizer)
		end
	end
	@require GLPK = "60bf3e95-4087-53dc-ae20-288a0d20c6a6" begin
		function add_GLPK!(m)
			JuMP.set_optimizer(m, GLPK.Optimizer)
		end
	end
	@require Xpress = "9e70acf3-d6c9-5be6-b5bd-4e2c73e3e054" begin
		function add_Xpress!(m)
			JuMP.set_optimizer(m, Gurobi.Optimizer)
		end
	end
	@require KNITRO = "67920dd8-b58e-52a8-8622-53c4cffbe346" begin
		function add_KNITRO!(m)
			JuMP.set_optimizer(m, KNITRO.Optimizer)
		end
	end
	@require Juniper = "2ddba703-00a4-53a7-87a5-e8b9971dde84" begin
		function add_Juniper!(m)
			JuMP.set_optimizer(m, Juniper.Optimizer)
		end
	end
	@require SCIP = "bdc5fa1a-e2a0-4ec3-bc77-c59bd234be9e" begin
		function add_SCIP!(m)
			JuMP.set_optimizer(m, SCIP.Optimizer)
		end
	end

end


export mean, StartATA, LoadSettings!, AddFriendSets!, AddEnemySets!, AddConstr!, AddOverlaps!, AddExpScore!, GroupByFriendSet!, AddObjFun!, Assemble!, PrintResults!, LoadDesign, RunATA! #mycopy,  fast_sort!, myqle#, optimize

ATA
end
