function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--first_sim"
            help = "an option with an argument"
            arg_type = Int64
            required = true
        "--last_sim"
            help = "another option with an argument"
            arg_type = Int64
            required = true
        "--save_dir"
            help = "another option with an argument"
            default = "/mnt/parscratch/users/bi1ahd/sim/simCSh_allo_d8/"
        "--param_file"
            help = "another option with an argument"
            default = "scripts/param_comb_ct_S_h_d4.arrow"
        "--K"
            help = "Carrying capacity"
            #arg_type = Float64
            default = 10.0
        "--c"
            help = "Predator interference"
            default = 0.0
        "--h"
            help = "Hill exponent of the functional response"
            default = 2.0
        "--K_corrected"
            help = "an option without argument, i.e. a flag"
            default = true
        "--d"
            default = nothing
        "--d_allometric_set"
            default = (ap = .4, ai = .4, ae = .4)
        "--tmax"
            help = "Number of timesteps"
            arg_type = Int
            default = 2000
        "--check_disconnected"
            help = "Should disconnected species be killed or removed?"
            arg_type = Bool
            default = true
        "--rebuild_after_disconnected"
            help = "Should the model be restarted without disconnected?"
            arg_type = Bool
            default = false
        "--re_run"
            help = "Should the model run until no more extinctions?"
            arg_type = Bool
            default = true
    end

    return parse_args(s)
end

function print_argument()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end
