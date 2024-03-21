# Trapezoid temperature function (TPF)

Trapezoid temperature function with four cardinal temperatures reflecting a range of optimum temperatures for photosynthesis[^1].

$\ f(T_{ab}) = ( T_{day} - T_{min}) / ( T_{optmin} - T_{min})$

$\ f(T_{cd}) = ( T_{max} - T_{day}) / ( T_{max} - T_{optmax}) $ 

Where:

$\ f(T) = 0;  T_{day} <= T_{min}$  or  $\ T_{day} >= T_{max} $

$\ f(T) = 1;  T_{optmin} <= T_{day} <= T_{optmax} $


[^1]: four parameters


Code for TPF:

=== "Python"

    ``` python
    # ----------------------------------------------
    # Trapezoidal Temperature Function (TPF) 
    # ----------------------------------------------
    def calculate_TPF(Tday, Tmin, Toptmin, Toptmax, Tmax):
        tpf = 0
        if (Toptmin > Toptmax):
            print("Min Optimum Temperature greater than Max Opt. Temperature")
            tpf = np.nan
        elif (Toptmax > Tmax):
            print("Max Optimum Temperature greater than Maximum Temperature")
            tpf = np.nan
        else:
            if ((Tday < Tmin) or (Tday > Tmax)):
                tpf = 0
            elif ((Tday >= Toptmin) and (Tday <= Toptmax)):
                tpf = 1
            elif (Tday < Toptmin):
                gradient, intercept, r_value, p_value, std_err = stats.linregress([Tmin,Toptmin],[0,1])
                tpf = Tday * gradient
            elif (Tday > Toptmax):
                gradient, intercept, r_value, p_value, std_err = stats.linregress([Toptmax, Tmax],[1,0])
                tpf = 1-((Tday-Toptmax)*abs(gradient))
        #
        return tpf
    ```

=== "R"

    ``` R
    # Trapezoidal Temperature Function
    calculateTPF <- function(Tday, Tmin, Toptmin, Toptmax, Tmax){
        tpf <- 0
        if (Toptmin > Toptmax){
            warning("Min Optimum Temperature greater than Maximum Opt. Temperature")
            tpf <- NA
        } else if (Toptmax > Tmax){
            warning("Max Optimum Temperature greater than Maximum Temperature")
            tpf <- NA
        } else {
            if ((Tday <= Tmin) | (Tday >= Tmax)){
            tpf <- 0
            } else if ((Tday >= Tmin) & (Tday < Toptmin)){
            interceptLeft <- lsfit(x=c(Tmin,Toptmin), y=c(0,1))$coefficients[1]
            slopeLeft <- lsfit(x=c(Tmin,Toptmin), y=c(0,1))$coefficients[2]
            tpf <- Tday * slopeLeft + interceptLeft
            } else if ((Tday > Toptmax) & (Tday <= Tmax)){
            interceptRight <- lsfit(x=c(Toptmax, Tmax), y=c(1,0))$coefficients[1]
            slopeRight <- lsfit(x=c(Toptmax, Tmax), y=c(1,0))$coefficients[2]
            tpf <- Tday * slopeRight + interceptRight
            } else if ((Tday >= Toptmin) & (Tday <= Toptmax)) { 
            tpf <- 1
            }
        }
        return(tpf)
    }
    calculateTPF_wrap <- function( Tday, Tmin, Toptmin, Toptmax, Tmax) {
        y <- vector("numeric", length (Tday)) 
        for (i in 1:length(Tday)) {
            y[[i]] <- calculateTPF( Tday=Tday[[i]], Tmin, Toptmin, Toptmax, Tmax) 
        }
        return(y)
    }
    ```
    The whole R code can be found in the Github repository at `Rcode` folder.

