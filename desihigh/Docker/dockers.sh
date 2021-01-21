# Create an image from the same dir as the Dockerfile. 
# docker build --tag desihigh:1.0.0 .

# Create a container from that image. 
# docker create desihigh:1.0.0
# docker create -t -i  572580b57b12 bash

# docker container ls --all

# -v /Users/MJWilson/Work/desi/legacy/docker/mjwilson:/src/mjwilson
# docker container run -it -p 8888:8888 74252bb1b6a1 /bin/bash

# Run in said container.                                                                                                                                                                                                                                                                                                     
# docker exec -ti c3a5b2fc7ec2 /bin/bash

# docker stop 25442f469369 
