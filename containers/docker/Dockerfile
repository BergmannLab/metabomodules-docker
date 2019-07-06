# Copyright 2019 Bergmann's Lab UNIL <mattia.tomasoni@unil.ch>
#
#    This file is part of METABOMODULES Tool.
#
#    METABOMODULES Tool is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    METABOMODULES Tool is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    GNU General Public License: <https://www.gnu.org/licenses/>.

FROM python:3
LABEL maintainer="Mattia Tomasoni <mattia.tomasoni@unil.ch>"

###############################################################################
# install octave and necessary python libs (python3 itself already installed)

RUN apt-get update && \
 apt-get -y install software-properties-common && \
 apt-get update && \
 apt-get -y install octave

RUN pip3 install numpy
RUN pip3 install sklearn
RUN pip3 install pandas
RUN pip3 install matplotlib

###############################################################################
# run method

ENTRYPOINT ["/bin/bash", "/metabomodules_code/metabomodules.sh"]
