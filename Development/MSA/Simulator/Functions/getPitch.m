function angle = getPitch(pitch_params, altitude)

    N = 100;
    initial_angle = pitch_params.first_angle;
    order = pitch_params.order;
    initial_altitude = pitch_params.initial_altitude;
    final_altitude = pitch_params.final_altitude;

    y = @(x) (-1)^order * initial_angle*(x).^order;

    x_vect = linspace(-1,0,N);
    y_vect = y(x_vect);

    new_vect = linspace(initial_altitude, final_altitude, N);

    if altitude < 11e3
        altitude = 11e3;
    end

    idx = sum(altitude >= new_vect);
    angle = y_vect(idx);
    
end

