FROM tensorflow/tensorflow:latest-gpu
MAINTAINER Matt Claudio (mc824@njit.edu)
# ENV NVIDIA_VISIBLE_DEVICES all
RUN pip3 install Keras sklearn matplotlib pandas yahoo_fin requests_html tensorflow Bio





