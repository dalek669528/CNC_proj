function fitness = checkIDfitness(data,sys)
    sim_out = sim(sys,data.Inputdata);
    fitness = 1-norm(data.Outputdata-sim_out)/norm(data.Outputdata-mean(data.Outputdata));
end