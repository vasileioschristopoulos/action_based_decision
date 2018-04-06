function converted_field=convertTwoDReferenceFramePolar(orig_field, fieldSize, NexpRewNeurons, effector_position, target_positions)

converted_field=zeros(1,fieldSize);
for j=1:size(target_positions,1)
    rx=target_positions(j,1);
    ry=target_positions(j,2);
    x_idx=round((rx*(90.0/NexpRewNeurons))+50.0)+1;
    y_idx=round(ry*2.0)+1;
    egocentric_pos=[rx-effector_position(1),ry-effector_position(2)];
    [effector_theta,effector_rho] = cart2pol(egocentric_pos(1),egocentric_pos(2));
    converted_field=converted_field+orig_field(x_idx,y_idx)*gauss(1:fieldSize, (effector_theta/pi)*180.0, 5);
end

