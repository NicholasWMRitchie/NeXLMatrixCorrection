using NeXLCore
using NeXLMatrixCorrection

k240 = NeXLCore.material("K240",Dict(n"O"=>0.340023, n"Mg"=>0.030154, n"Si"=>0.186986, n"Ti"=>0.059950, n"Zn"=>0.040168, n"Zr"=>0.074030, n"Ba"=>0.268689),missing)
sio2 = atomicfraction("Quartz",Dict(n"Si"=>1,n"O"=>2))
mgo = atomicfraction("MgO",Dict(n"Mg"=>1,n"O"=>1))
baf2 = atomicfraction("Barium Fluoride",Dict(n"Ba"=>1,n"F"=>2))
ti, zn, zr = pure(n"Ti"), pure(n"Zn"), pure(n"Zr")

e0, θ = 20.0, 40.0

zafK240Si = xppZAF(k240, n"Si K", e0, θ)
zafK240Mg = xppZAF(k240, n"Mg K", e0, θ)
zafK240Ba = xppZAF(k240, n"Ba L3", e0, θ)
zafK240Ti = xppZAF(k240, n"Ti K", e0, θ)
zafK240Zn = xppZAF(k240, n"Zn K", e0, θ)
zafK240Zr = xppZAF(k240, n"Zr K", e0, θ)
zafK240O = xppZAF(k240, n"O K", e0, θ)

zafSi = xppZAF(sio2, n"Si K", e0, θ)
zafMg = xppZAF(mgo, n"Mg K", e0, θ)
zafBa = xppZAF(baf2, n"Ba L3", e0, θ)
zafTi = xppZAF(ti, n"Ti K", e0, θ)
zafZn = xppZAF(zn, n"Zn K", e0, θ)
zafZr = xppZAF(zr, n"Zr K", e0, θ)
zafO = xppZAF(sio2, n"O K", e0, θ)

summarize(Dict(zafK240Si=>zafSi, zafK240Mg=>zafMg, zafK240Ba=>zafBa, zafK240Ti=>zafTi, zafK240Zn=>zafZn, zafK240Zr=>zafZr, zafK240O=>zafO))
