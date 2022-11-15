using BEFWM2
using Plots
using DataFrames
using CSV

A = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0] 
foodweb = FoodWeb(A)

# Functional response
## Preference of consumer 
myω = zeros(4, 4)
myω[:,1] = [0, 1, .98, 0]
myω[4,:] = [0, .92, 1 - .92, 0]
## Predator interference
myc = repeat([0], 4)
myc

# Half-saturation constant
myB0 = [0, 0.16129, .90, .5]

# Biological rates
## Assimilitation efficiency or ingestion rates (in Vasseur & Fox, 2007)
## Efficiency is implicitly one for all consumers in Vasseur & Fox
mye = zeros(4, 4)
# Consumer
mye[:,1] = [0, 1, 1, 0] 
# Predator
mye[4,:] = [0, 1, 1, 0]

# Metabolic rates and maximum ingestion rates 
## From McCann (1998):
x = [0, .40, .20, .08]
y = [0, 2.009, 3.50, 5.0]

#  

####################
#  McCann version  #
####################
# Define x and y

bioener = BioenergeticResponse(foodweb,
                               h = 1,
                               # Half saturation-constant
                               B0 = myB0,
                               # Consumer preference
                               ω = myω,
                               # Predator interference
                               c = myc
                              )

biorate = BioRates(foodweb,
        r = [1.0, 0, 0, 0],
        e = mye,
        x = x, 
        y = y
       )

params = ModelParameters(foodweb,
                functional_response = bioener,
                biorates = biorate,
                environment = Environment(foodweb, K = 1.0)
               )

m = simulate(params, rand(4))
df = DataFrame(m)
CSV.write("./mcann_original.csv",df)

mccann_fig = plot(m)


###########################
#  McCann version with Z  #
###########################
# myω

fb = FoodWeb(A, Z = 9.68)
bner = BioenergeticResponse(fb,
                               h = 1,
                               # Half saturation-constant
                               B0 = myB0,
                               # Consumer preference
                               ω = myω,
                               # Predator interference
                               c = myc
                              )

br = BioRates(fb, r = [1.0, 0, 0, 0], e = mye)

pb = ModelParameters(fb,
                functional_response = bner,
                biorates = br,
                environment = Environment(fb, K = 1.0)
               )

mccann_z_sim = simulate(pb, rand(4))
df = DataFrame(mccann_z_sim)
CSV.write("mcann_z.csv",df)

mccann_z_fig = plot(mccann_z_sim)
mccann_z_fig

###########
#  Brose  #
###########
# Do not define x and y but Z 
# consumer preference even
#
foodwebb = FoodWeb(A, Z = 4.94)
bioenerb = BioenergeticResponse(foodwebb,
                               h = 1,
                               # Half saturation-constant
                               B0 = [0, .5, .5, .5],
                              )
bioenerb.ω

biorateb = BioRates(foodwebb, r = [1.0, 0, 0, 0], e = mye)
biorate

paramsb = ModelParameters(foodwebb,
                functional_response = bioenerb,
                biorates = biorateb,
                environment = Environment(foodwebb, K = 1.0)
               )

mb = simulate(paramsb, rand(4))
df = DataFrame(mb)
CSV.write("brose_little_z.csv",df)

brose_fig = plot(mb)
brose_fig

#################
#  BROSE big Z  #
#################
#Z corresponding to the weak interaction P-C2 in McCann 

foodwebb = FoodWeb(A, Z = 242)
bioenerb = BioenergeticResponse(foodwebb,
                               h = 1,
                               # Half saturation-constant
                               B0 = [0, 0.5, .5, .5],
                              )

biorateb = BioRates(foodwebb, r = [1.0, 0, 0, 0], e = mye)
biorate

paramsb = ModelParameters(foodwebb,
                functional_response = bioenerb,
                biorates = biorateb,
                environment = Environment(foodwebb, K = 1.0)
               )

mb = simulate(paramsb, rand(4))
df = DataFrame(mb)
CSV.write("brose_big_z.csv",df)

brose_fig_Z_weak = plot(mb)




