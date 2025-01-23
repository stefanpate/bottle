# jre-build
FROM eclipse-temurin:17 AS jre-build

# Create a custom Java runtime
RUN $JAVA_HOME/bin/jlink \
         --add-modules java.se \
         --strip-debug \
         --no-man-pages \
         --no-header-files \
         --compress=2 \
         --output /javaruntime

# chemaxon jchem (includes cxcalc)
COPY ./jchem_linux_24.3.2.deb jchem.deb
RUN apt-get install ./jchem.deb

# --- runtime layer
FROM python:3.12-slim-bookworm AS runtime

ENV JAVA_HOME=/opt/java/openjdk
ENV CHEMAXON_HOME=/opt/chemaxon/jchemsuite
ENV PATH="${JAVA_HOME}/bin:${CHEMAXON_HOME}/bin:${PATH}"

COPY --from=jre-build /javaruntime ${JAVA_HOME}
COPY --from=jre-build /opt/chemaxon/jchemsuite ${CHEMAXON_HOME}
