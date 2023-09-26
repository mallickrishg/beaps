function [ux,uz] = compute_bemdisplacements(rcv,obs,slip_d,slip_n)

[Gdx,Gdz,Gnx,Gnz] = geometry.computeDisplacementKernels(rcv,obs);
ux = Gdx*slip_d + Gnx*slip_n;
uz = Gdz*slip_d + Gnz*slip_n;


end