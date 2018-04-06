function direction = atan2d(y,x)
% Function calculates the arc tangent of y/x and places the result in
% range of [0..360]
direction = 0;
if x == 0
    if y == 0
    elseif y > 0
        direction = 90;
    else
        direction = 270;
    end
elseif x > 0
    if y >= 0
        direction = atand(y/x);
    else
        direction = atand(y/x) + 360;
    end
elseif x < 0
    if y == 0
        direction = 180;
    else
        direction = atand(y/x) + 180;
    end
end
end

