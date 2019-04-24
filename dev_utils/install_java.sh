#!/bin/bash

sudo mkdir ~/Downloads/
cd ~/Downloads/
sudo wget http://storage.googleapis.com/travis-java/jdk-8u211-linux-i586.tar.gz
if [ -d "/usr/lib/jvm" ]; then rm -rf /usr/lib/jvm; fi
sudo mkdir /usr/lib/jvm
cd /usr/lib/jvm
sudo tar -xvzf ~/Downloads/jdk-8u211-linux-i586.tar.gz

echo "PATH=\"/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/usr/lib/jvm/jdk1.8.0_211/bin:/usr/lib/jvm/jdk1.8.0_211/db/bin:/usr/lib/jvm/jdk1.8.0_211/jre/bin\"
J2SDKDIR=\"/usr/lib/jvm/jdk1.8.0_211\"
J2REDIR=\"/usr/lib/jvm/jdk1.8.0_211/jre*\"
JAVA_HOME=\"/usr/lib/jvm/jdk1.8.0_211\"
DERBY_HOME=\"/usr/lib/jvm/jdk1.8.0_211/db\"" > /etc/environment

source /etc/environment

sudo update-alternatives --install "/usr/bin/java" "java" "/usr/lib/jvm/jdk1.8.0_211/bin/java" 0
sudo update-alternatives --install "/usr/bin/javac" "javac" "/usr/lib/jvm/jdk1.8.0_211/bin/javac" 0
sudo update-alternatives --set java /usr/lib/jvm/jdk1.8.0_211/bin/java
sudo update-alternatives --set javac /usr/lib/jvm/jdk1.8.0_211/bin/javac
