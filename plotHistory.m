function plotHistory(history, end_time, target_positions, cue_color, event_times)

event1_time=50;
event2_time=100;

for i=1:length(history.record_trials)

    movement_time=min([min(find(history.reaching_movements(i,:,2)>0)) min(find(history.reaching_movements(i,:,2)<0)) min(find(history.reaching_movements(i,:,1)>0)) min(find(history.reaching_movements(i,:,2)<0)) min(find(history.saccade_movements(i,:,2)>0)) min(find(history.saccade_movements(i,:,2)<0)) min(find(history.saccade_movements(i,:,1)>0)) min(find(history.saccade_movements(i,:,2)<0))]);

    figure(2+(i-1)*11+1);    
    subplot(7, 1, 1);
    plot2dField(history.dnf_stim_eye(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus eye',0);
    subplot(7, 1, 2);
    plot2dField(history.dnf_stim_hand(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus hand',0);
    subplot(7, 1, 3);
    plot2dField(history.cue(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[0,1.2],'cue',0);
    subplot(7, 1, 4);
    plot2dField(history.effort_sac(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-.2,2],'saccade effort',0);
    subplot(7, 1, 5);
    plot2dField(history.effort_rch(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-.2,2],'reach effort',0);
    subplot(7, 1, 6);
    plot2dField(history.dnf_sac(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'saccade u',0);
    subplot(7, 1, 7);
    plot2dField(history.dnf_rch(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'reach u',0);
    
    field=figure(2+(i-1)*11+2);
    %subplot(2, 2, 1);
    stim_eye=squeeze(history.dnf_stim_eye(i,1:end_time,:))';
    stim_eye(find(stim_eye<-10))=-10;
    plot3dField(stim_eye,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus eye',0);
    saveas(field,['trial_' num2str(i) '.field.stim_eye.ai'], 'ai');
    field=figure(2+(i-1)*11+3);
    %subplot(2, 2, 2);
    stim_hand=squeeze(history.dnf_stim_hand(i,1:end_time,:))';
    stim_hand(find(stim_hand<-10))=-10;
    plot3dField(stim_hand,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus hand',0);
    saveas(field,['trial_' num2str(i) '.field.stim_hand.ai'], 'ai');
    field=figure(2+(i-1)*11+4);
    %subplot(2, 2, 3);
    dnf_eye=squeeze(history.dnf_sac(i,1:end_time,:))';
    dnf_eye(find(dnf_eye<-10))=-10;
    plot3dField(dnf_eye,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'saccade u',0);
    saveas(field,['trial_' num2str(i) '.field.dnf_eye.ai'], 'ai');
    field=figure(2+(i-1)*11+5);
    %subplot(2, 2, 4);
    dnf_hand=squeeze(history.dnf_rch(i,1:end_time,:))';
    dnf_hand(find(dnf_hand<-10))=-10;
    plot3dField(dnf_hand,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'reach u',0);
    saveas(field,['trial_' num2str(i) '.field.dnf_rch.ai'], 'ai');
    %plot2svg(['trial_' num2str(i) '.field.svg'], field, 'png');
    %saveas(field,['trial_' num2str(i) '.field.eps'], 'eps');
    %print('-dpdf','-painters', ['trial_' num2str(i) '.field.pdf']);
    
    figure(2+(i-1)*11+6);
    subplot(7, 1, 1);
    plot2dField(history.dnf_stim_eye(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus eye',1);
    subplot(7, 1, 2);
    plot2dField(history.dnf_stim_hand(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus hand',1);
    subplot(7, 1, 3);
    plot2dField(history.cue(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[0,1.2],'cue',1);
    subplot(7, 1, 4);
    plot2dField(history.effort_sac(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-.2,2],'saccade effort',1);
    subplot(7, 1, 5);
    plot2dField(history.effort_rch(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-.2,2],'reach effort',1);
    subplot(7, 1, 6);
    plot2dField(history.dnf_sac(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'saccade u',1);
    subplot(7, 1, 7);
    plot2dField(history.dnf_rch(i,1:end_time,:),history.record_trials(i),event1_time,event2_time,movement_time,event_times,[-10,10],'reach u',1);

    field_gray=figure(2+(i-1)*11+7);
    %subplot(2, 2, 1);
    plot3dField(stim_eye,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus eye',1);
    saveas(field_gray,['trial_' num2str(i) '.field_gray.stim_eye.ai'], 'ai');
    %subplot(2, 2, 2);
    field_gray=figure(2+(i-1)*11+8);
    plot3dField(stim_hand,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'stimulus hand',1);
    saveas(field_gray,['trial_' num2str(i) '.field_gray.stim_hand.ai'], 'ai');
    %subplot(2, 2, 3);
    field_gray=figure(2+(i-1)*11+9);
    plot3dField(dnf_eye,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'saccade u',1);
    saveas(field_gray,['trial_' num2str(i) '.field_gray.dnf_eye.ai'], 'ai');
    %subplot(2, 2, 4);
    field_gray=figure(2+(i-1)*11+10);
    plot3dField(dnf_hand,history.record_trials(i),end_time,event1_time,event2_time,movement_time,event_times,[-10,10],'reach u',1);
    saveas(field_gray,['trial_' num2str(i) '.field_gray.dnf_hand.ai'], 'ai');
    %plot2svg(['trial_' num2str(i) '.field_gray.svg'], field_gray, 'png');
    %saveas(field_gray,['trial_' num2str(i) '.field_gray.eps'], 'eps');
    %print('-depsc','-painters', ['trial_' num2str(i) '.field_gray.eps']);

    traj=figure(2+(i-1)*11+11);
    plot(0,0,'.');
    hold on;
    symbol='.';
    if history.dopamine(i)>0
        symbol='o';
    end
    plot(squeeze(history.reaching_movements(i,:,1)), squeeze(history.reaching_movements(i,:,2)), symbol, 'color', 'g');
    plot(squeeze(history.saccade_movements(i,:,1)), squeeze(history.saccade_movements(i,:,2)), symbol, 'color', 'r');
    for j=1:size(target_positions,1)
        rx=target_positions(j,1);
        ry=target_positions(j,2);
        plot(rx,ry,'ko','MarkerSize',10,'MarkerFaceColor','k');
        myCircle(5 ,rx,ry,100,'k');
    end
    if length(cue_color)>0
        plot(0,21,[cue_color 'o'],'MarkerSize',10,'MarkerFaceColor',cue_color);
    end
    for j=1:length(event_times)
        rx=0;
        ry=0;
        if history.reaching_movements(i,event_times(j),1)~=0 || history.reaching_movements(i,event_times(j),2)~=0
            rx=history.reaching_movements(i,event_times(j),1);
            ry=history.reaching_movements(i,event_times(j),2);
        elseif history.saccade_movements(i,event_times(j),1)~=0 || history.saccade_movements(i,event_times(j),2)~=0
            rx=history.saccade_movements(i,event_times(j),1);
            ry=history.saccade_movements(i,event_times(j),2);
        end
        plot(rx,ry,'ro--','MarkerSize',15);
    end
    hold off;
    set(gca,'ylim',[0,50],'xlim',[-45,45]);
    daspect([1 1 1]);
    saveas(traj,['trial_' num2str(i) '.trajectory.ai'], 'ai');
end

max_cue_w=max([max(history.cue_sac_w(:)) max(history.cue_rch_w(:))]);
min_cue_w=min([min(history.cue_sac_w(:)) min(history.cue_rch_w(:))]);
figure(2+length(history.record_trials)*11+1);
subplot(5,1,1);
imagesc(squeeze(history.cue_sac_w(:,1,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 1 -> saccade');
colorbar();
subplot(5,1,2);
imagesc(squeeze(history.cue_sac_w(:,2,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 2 -> saccade');
colorbar();
subplot(5,1,3);
imagesc(squeeze(history.cue_rch_w(:,1,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 1 -> reach');
colorbar();
subplot(5,1,4);
imagesc(squeeze(history.cue_rch_w(:,2,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 2 -> reach');
colorbar();
subplot(5,1,5);
%imagesc(history.expected_rew_w');
%xlabel('trial');
%ylabel('expected reward');
%subplot(4,1,4);
plot(history.dopamine);
xlabel('trial');
ylabel('dopamine');
ylim([0 1.1]);

figure(2+length(history.record_trials)*11+2);
subplot(5,1,1);
imagesc(squeeze(history.cue_sac_w(:,1,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 1 -> saccade');
colormap(gray);
colorbar();
subplot(5,1,2);
imagesc(squeeze(history.cue_sac_w(:,2,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 2 -> saccade');
colormap(gray);
colorbar();
subplot(5,1,3);
imagesc(squeeze(history.cue_rch_w(:,1,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 1 -> reach');
colormap(gray);
colorbar();
subplot(5,1,4);
imagesc(squeeze(history.cue_rch_w(:,2,:))',[min_cue_w,max_cue_w]);
xlabel('trial');
ylabel('cue 2 -> reach');
colormap(gray);
colorbar();
subplot(5,1,5);
%imagesc(history.expected_rew_w');
%xlabel('trial');
%ylabel('expected reward');
%subplot(4,1,4);
plot(history.dopamine);
xlabel('trial');
ylabel('dopamine');
ylim([0 1.1]);

er_w=zeros(size(history.expected_rew_w,1),size(history.dnf_stim_eye,3));
for i=1:size(history.expected_rew_w,1)
    er_w(i,:)=convertTwoDReferenceFramePolar(squeeze(history.expected_rew_w(i,:,:)),size(history.dnf_stim_eye,3),size(history.expected_rew_w,2),[0 0],target_positions);
end
figure(2+length(history.record_trials)*11+3);
subplot(2,1,1);
imagesc(er_w');
xlabel('trial');
ylabel('expected reward');
colorbar();
subplot(2,1,2);
plot(history.dopamine);
xlabel('trial');
ylabel('dopamine');
ylim([0 1.1]);

figure(2+length(history.record_trials)*11+4);
subplot(2,1,1);
imagesc(er_w');
xlabel('trial');
ylabel('expected reward');
colormap(gray);
colorbar();
subplot(2,1,2);
plot(history.dopamine);
xlabel('trial');
ylabel('dopamine');
ylim([0 1.1]);
end

function plot2dField(field,trial_num,event1_time,event2_time,movement_time,event_times,range,y_axis,gray)
    c='r';
    if gray>0
        c='w';
    end
    title(['Trial ' num2str(trial_num)]);
    imagesc(squeeze(field)',range);
    hold on
    plot([event1_time,event1_time],[1,181],[c '--']);
    plot([event2_time,event2_time],[1,181],[c '--']);
    plot([movement_time,movement_time],[1,181],[c '--']);
    for j=1:length(event_times)
        plot([event_times(j),event_times(j)],[1,181],[c '--']);
    end
    hold off
    xlabel('time');
    ylabel(y_axis);
    if gray>0
        colormap('gray');
    end
    colorbar();
end

function plot3dField(field,trial_num,end_time,event1_time,event2_time,movement_time,event_times,range,y_axis,gray)
    c='r';
    if gray>0
        c='w';
    end
    title(['Trial ' num2str(trial_num)]);
    surf([1:end_time],[1:181],field, 'EdgeColor', 'none');
    hold on;
    plot3(ones(1,181)*(event1_time-1),[1:181],field(:,(event1_time-1))+.5,[c '--']);
    plot3(ones(1,181)*(event2_time-1),[1:181],field(:,(event2_time-1))+.5,[c '--']);
    plot3(ones(1,181)*(movement_time-1),[1:181],field(:,(movement_time-1))+.5,[c '--']);
    for j=1:length(event_times)
        plot3(ones(1,181)*(event_times(j)-1),[1:181],field(:,(event_times(j)-1))+.5,[c '--']);
    end
    view([-51.5,42]);
    ylim([1,181]);
    zlim(range);
    caxis(range);
    xlabel('time');
    ylabel(y_axis);
    if gray>0
        colormap('gray');
    end
    colorbar();
end
    
