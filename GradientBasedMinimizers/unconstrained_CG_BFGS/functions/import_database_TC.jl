#= ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin Julia
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#

"""
    find_end(start, lines)

    Finds the first empty line starting "start"

"""
 function find_end(start, lines; case = 0)
    ending  = -1
    i       = start
    success = 0
    while success != 1
        i += 1
        
        if length(strip(lines[i])) == 0
            success = 1
        end

        if case == 1
            tmp =  match.(r"\w(?=\()", lines[i+1])
            if ~isnothing(tmp)
                if strip(tmp.match) == "v"
                    success = 0
                end
            end
        end


    end
    ending = i

    return ending
end


"""
    import_activity_models_TC(filename::String)

    This function reads the activity model from the NewForm output of THERMOCALC (for MAGEMin)

"""
function import_activity_models_TC(filename::String)

    lines   = readlines(filename)
    n_lines = length(lines)

    tmp_ss  = Int64[]
    start   = 1

    # here we look for the begining of the solution phase model definitions
    for i = 2:n_lines
        if lines[i] == "#" && lines[i-1] == "#"
            start = i
            break;
        end
    end

    # here we look for the begining of the definition for every solution phase model
    for i = start:n_lines
        if length(lines[i]) >= 3 && length(lines[i-1]) >= 1
            if lines[i][1:3] == " ==" && lines[i-1][:] == "#"
                push!(tmp_ss,i)
            end
        end
    end


    ss = Dict()
    for i = 1:length(tmp_ss)-1

        sg, sf, ni, im, mk, pr, la, sfcv = [], [], [], [], [], [], [], [] 
        sg_end, sf_end, ni_end, im_end, mk_end, pr_end, la_end, sfcv_end = [], [], [], [], [], [], [], [] 
        ss_name = String

        for j = tmp_ss[i]:tmp_ss[i+1]
            if lines[j] == " starting guesses"
                sg = j
                sg_end = find_end(sg, lines)

                ss_name = String(match.(r"\((.*?)\)", lines[j+1])[1])
                ss[ss_name] = Dict( "labels"    => [], #name, abr, var
                                    "sf"        => [],
                                    "cv"        => [],
                                    "p"         => [],
                                    "bounds"    => [], #lower bounds, upper bounds
                                    "sf2cv"     => [],
                                    "idm"       => [],
                                    "W"         => [],
                                    "v"         => [], #if v, then symmetry = 0
                                    "make"      => [],
                                       )
            end

            if lines[j] == " sf2cv"
                sfcv = j
                sfcv_end = find_end(sfcv, lines)
            end
            if lines[j] == " site fractions"
                sf = j
                sf_end = find_end(sf, lines)
            end
            if lines[j] == " non-ideality by symmetric formalism" || lines[j] == " non-ideality by van laar"
                ni = j
                ni_end = find_end(ni, lines; case = 1)
            end
            if lines[j] == " ideal mixing activities"
                im = j
                im_end = find_end(im, lines)
            end
            if lines[j] == " \"make\" end-members"
                mk = j
                mk_end = find_end(mk, lines)
            end
            if lines[j] == " proportions"
                pr = j
                pr_end = find_end(pr, lines)
            end
            if lines[j] == " labels"
                la = j
                la_end = find_end(la, lines)
            end

        end

        # retrieve the margules W's parameters:
        # print(ni," ",ni_end,"\n")
        for j = ni+1:ni_end
            if length(lines[j]) > 0

                var = match.(r"\w(?=\()", lines[j]).match
                if var == "W"
                    fields = split(lines[j],"=")

                    if length(fields) == 2
                        tmp = (strip(fields[1]),strip(fields[2]))
                        push!(ss[ss_name]["W"],tmp)
                    else
                        print("something went wrong for $ss_name, splitting the line for Ws should results in 2 terms: t1 = t2\n")
                    end
                end
                if var == "v"
                    fields = split(lines[j],"=")

                    if length(fields) == 2
                        tmp = (strip(fields[1]),strip(fields[2]))
                        push!(ss[ss_name]["v"],tmp)
                    else
                        print("something went wrong for $ss_name, splitting the line for vs should results in 2 terms: t1 = t2\n")
                    end
                end

            end
        end
        
        # retrieve the labeling
        if isempty(la)
            print("Labels are not provided for $ss_name or formatting is wrong\n")
        else
            # n_lbl = (sf-1)-(la+1)
            for j = la+1:la_end
                if length(lines[j]) > 0
                    fields = split(lines[j],":")
                    if length(fields) == 3
                        tmp = (strip(fields[1]),strip(fields[2]),strip(fields[3]))
                        
                        push!(ss[ss_name]["labels"],tmp)
                    else
                        print("something went wrong for $ss_name, splitting the line for labels should results in 3 terms: t1 : t2 : t3\n")
                    end
                end
            end

        end

        if ~isempty(sfcv)
            # retrieve site fractions expressions
            for j = sfcv+1:sfcv_end
                if length(lines[j]) > 0
                    fields = split(lines[j],"->")

                    if length(fields) == 2
                        tmp = (strip(fields[1]),strip(fields[2]))
                        push!(ss[ss_name]["sf2cv"],tmp)
                    else
                        print("something went wrong for $ss_name, splitting the line for sf2cv should results in 2 terms: t1 = t2\n")
                    end
                end
            end
        end


        # retrieve site fractions expressions
        for j = sf+1:sf_end
            if length(lines[j]) > 0
                fields = split(lines[j],"=")

                if length(fields) == 2
                    tmp = (strip(fields[1]),strip(fields[2]))
                    push!(ss[ss_name]["sf"],tmp)
                else
                    print("something went wrong for $ss_name, splitting the line for sf should results in 2 terms: t1 = t2\n")
                end
            end
        end

        # retrieve compositional variables
        for j = sg+1:sg_end
            if length(lines[j]) > 0
                fields = split(lines[j],"(")

                if length(fields) == 2
                    tmp = strip(fields[1])
                    push!(ss[ss_name]["cv"],tmp)
                    fields_b = split(fields[2],"range")
                    fields_c = split(fields_b[2],"<>")
                    tmp = (strip(fields_c[1]),strip(fields_c[2]))
                    push!(ss[ss_name]["bounds"],tmp)
                else
                    print("something went wrong for $ss_name, splitting the line for cv should results in 2 terms: t1 ( t2\n")
                end
            end
        end

        # retrieve end-member fractions expressions
        for j = pr+1:pr_end
            if length(lines[j]) > 0
                fields = split(lines[j],"=")

                if length(fields) == 2
                    tmp = (strip(fields[1]),strip(fields[2]))
                    push!(ss[ss_name]["p"],tmp)
                else
                    print("something went wrong for $ss_name, splitting the line for end-member fractions should results in 2 terms: t1 = t2\n")
                end
            end
        end

        # retrieve idm expressions
        for j = im+1:im_end
            if length(lines[j]) > 0
                fields = split(lines[j],"=")

                if length(fields) == 2
                    tmp = (strip(fields[1]),strip(fields[2]))
                    tmp2 = replace.(tmp, r"\*\*" => "^")
                    push!(ss[ss_name]["idm"],tmp2)
                else
                    print("something went wrong for $ss_name, splitting the line for idm should results in 2 terms: t1 = t2\n")
                end
            end
        end

        # retrieve endmembers makes
        if isempty(mk)
            for j=1:length(ss[ss_name]["p"])
                tmp = (ss[ss_name]["p"][j][1],ss[ss_name]["p"][j][1])
                push!(ss[ss_name]["make"],tmp)
            end

        else  # case when there is composite endmember (we overide)
            for j=1:length(ss[ss_name]["p"])
                lhs = ss[ss_name]["p"][j][1]
                tmp = (ss[ss_name]["p"][j][1],ss[ss_name]["p"][j][1])
                for k = mk+1:mk_end
                    
                    if length(lines[k]) > 0
                        fields  = split(lines[k],"=")
                        make    = strip(fields[1])
                        if lhs == make
                            rhs     = strip(replace.( fields[2], r"\([^)]*\)"           => "")      )
                            rhs2    = replace(              rhs, r"\b(\w+)-(\w+)\b"     => s"\1_\2" )
                            rhs3    = replace(             rhs2, r"(\d+)\s+([a-zA-Z]+)" => s"\1*\2" )

                            tmp = (ss[ss_name]["p"][j][1],rhs3)
                        end
                    end  
                end
                push!(ss[ss_name]["make"],tmp)
            end
        end

    end


    return ss
end

ss = import_activity_models_TC("TC_database/igneous_set_full_descriptions_NewForm.txt")
