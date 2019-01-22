#!/bin/bash


if sudo apt-get install oracle-java8-installer -y ; then
    echo "Java8 installed."
else 
    url='http://javadl-esd-secure.oracle.com/update/baseline.version'
    version=$(curl -L $url | grep 1.8)

    j_dir="jdk$version"
    v_num=$(echo $version | cut -d "_" -f 2)
    download_link=$(curl -i https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html | grep jdk-8u$v_num-linux-x64.tar.gz | cut -d '"' -f12 | cut -d '/' -f1-8)

    link=$(curl -i https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html | grep -o https://download.oracle.com/.*/jdk-8u$v_num-linux-x64.tar.gz)
    link=${link/https/http}

    cd /var/cache/oracle-jdk8-installer
    checksum=$(sha256sum jdk-8u201-linux-x64.tar.gz | cut -d ' ' -f 1)

    sed_string='s|JAVA_VERSION=8u.*|JAVA_VERSION=8u'$v_num'|'
    cd /var/lib/dpkg/info && sudo sed -i $sed_string oracle-java8-installer.*

    sed_string='s|PARTNER_URL=http://download.oracle.com.*|PARTNER_URL='$download_link'/|'
    echo $sed_string
    sudo sed -i $sed_string oracle-java8-installer.*

    sed_string='s|SHA256SUM_TGZ=.*|SHA256SUM_TGZ='\"$checksum\"'|'
    sudo sed -i $sed_string oracle-java8-installer.*

    sed_string='s|J_DIR=jdk.*|J_DIR='$j_dir'|'
    sudo sed -i $sed_string oracle-java8-installer.* 
fi
