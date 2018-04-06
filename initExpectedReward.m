function expected_reward=initExpectedReward(NexpRewNeurons,protocolParams)

expected_reward=zeros(NexpRewNeurons,NexpRewNeurons);
for i=1:size(protocolParams.target_positions,1)
    rx=protocolParams.target_positions(i,1);
    ry=protocolParams.target_positions(i,2);
    for x=1:NexpRewNeurons
        for y=1:NexpRewNeurons
            % x coordinate for index x, eligspatial is 100x100 and codes for x=-45 to x=45
            x_coord=(x-NexpRewNeurons/2.0+1.0)/(NexpRewNeurons)*(45.0--45.0);
            % y coordinate for index y, eligspatial is 100x100 and codes for y=0 to x=50
            y_coord=(y-1.0)/(NexpRewNeurons-1)*(50-0.0);
            expected_reward(x,y)=expected_reward(x,y)+protocolParams.init_expected_reward(i)*gauss(x_coord, rx, 5)*gauss(y_coord, ry, 5);
        end
    end
end
expected_reward = expected_reward/max(expected_reward(:));
