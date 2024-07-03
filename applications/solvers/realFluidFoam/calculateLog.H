/*
scalar vaporPenetration = 0.;

label fuel_index = 0;
if (composition.species().found("N2"))
{
    fuel_index = composition.species()["N2"];
}

vector position(1.5, 1, 1.7);
vector direction(0, 0, -1);

scalarField Y_fuel(T.size(), 0.);
Y_fuel = Y[fuel_index];

forAll (Y_fuel,cellI)
{
    if (Y_fuel[cellI] >= 0.001)
    {
        vector raw = position - mesh.C()[cellI];
        if (mag(raw&direction) > vaporPenetration)
        {
            vaporPenetration = mag(raw&direction);
        }
    }
}

reduce(vaporPenetration, maxOp<scalar>());


label OH_index = 0;
if (composition.species().found("OH"))
{
    OH_index = composition.species()["OH"];
}
scalar Y_OH_max = gMax(Y[OH_index]);

*/


//scalar T_mean = T.weightedAverage(mesh.V()).value();
scalar T_max = gMax(T);
scalar T_min = gMin(T);
scalar p_max = gMax(p);
scalar p_min = gMin(p);
//scalar T_max_old = gMax(T.oldTime());
//scalar dTdt = (T_max-T_max_old)/runTime.deltaTValue();


logSummaryFile
    << runTime.timeName() << tab		
//    << p.weightedAverage(mesh.V()).value() << tab
//    << T.weightedAverage(mesh.V()).value() << tab
//    << rho.weightedAverage(mesh.V()).value() << tab
//    << Qdot.weightedAverage(mesh.V()).value() << tab
//    << Qdot.average().value() << tab
//    << composition.Y("C12H26").weightedAverage(mesh.V()).value() << tab
    << rho.weightedAverage(mesh.V()).value() << tab	
    << p_min << tab
    << p_max << tab
    << T_min << tab
    << T_max << tab
//    << Y_OH_max << tab
//    << vaporPenetration << tab
//    << (sqrt((2.0/3.0)*turbulence->k()))().weightedAverage(mesh.V()).value() << tab
    << runTime.elapsedCpuTime() << tab
    << endl;

