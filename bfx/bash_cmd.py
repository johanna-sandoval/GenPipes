#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from core.job import Job

def mkdir(folder, remove=False):
    return Job(
        [],
        [folder],
        [],
        command="""\
mkdir -p {directory}""".format(
            directory=folder
        ),
        removable_files=[folder] if remove else []
    )

def chgdir(folder):
    return Job(
        [],
        [folder],
        [],
        command="""\
cd {directory}""".format(
            directory=folder
        ),
    )

def ln(target_file, link, out_dir=None):
    folder = os.path.dirname(link)
    return Job(
        [target_file],
        [link],
        [],
        command="""\
ln -s -f \\
  {target_file} \\
  {link}""".format(
            #target_file=os.path.relpath(target_file, folder),
            target_file=os.path.abspath(os.path.join(out_dir, target_file)),
            link=os.path.abspath(os.path.join(out_dir, link)),
            folder=folder
        ),
        removable_files=[link]
    )

def mv(source, target):
    return Job(
        [source],
        [target],
        [],
        command="""\
mv {source} \\
   {dest}""".format(
            source=source,
            dest=target
        )
    )

def rm(source):
    return Job(
        [source],
        [],
        command="""\
rm -rf {source}""".format(
            source=source,
        )
    )