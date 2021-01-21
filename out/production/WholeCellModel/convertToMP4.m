clc;
clear all;
close all;
% Browse Video File :
[ video_file_name,video_file_path ] = uigetfile({'/Users/shicheng/Desktop/invagination/id-61.mov'});
if(video_file_path==0)
    return;
end
% Output path
output_image_path = fullfile(video_file_path,[video_file_name(1:strfind(video_file_name,'.')-1),'.mp4']);
% mkdir(output_image_path);
input_video_file = [video_file_path,video_file_name];
% Read Video
videoFReader = VideoReader(input_video_file);
% Write Video
videoFWrite = VideoWriter(output_image_path,'MPEG-4');
open(videoFWrite);
final = abs(videoFReader.Duration*videoFReader.FrameRate)
if(final > 2e6)
    final = 2e6
end
for count = 1:250:final
    disp(count);
    key_frame = read(videoFReader,count);
    %key_frame = flip(key_frame,2);
    writeVideo(videoFWrite,key_frame);
end
% Release video object
% key_frame = readframe(videoFReader,final);
close(videoFReader);
close(videoFWrite);
disp('COMPLETED... (-_-)');