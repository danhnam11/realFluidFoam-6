//Info<< "Total cylinder mass: " << fvc::domainIntegrate(rho).value() << endl;

OFstream logSummaryFile
(
    runTime.path()/("logSummary." + runTime.timeName() + ".dat")
);

logSummaryFile
<< "variables=" << '\n'
<< "Time " << "rho(VolumeAveraged) " << "p_min(Pa) " << "p_max(Pa) "  << "T_min(K) " << "T_max(K) " << "ExecutionTime " <<  endl;
