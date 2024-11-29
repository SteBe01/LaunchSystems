function angle = getPitch(params, altitude)

    if altitude < 11e3
        altitude = 11e3;
    end

    data = params.angle_data;

    idx = sum(altitude >= data(:,1));
    angle = data(idx,2);
    
end

