docker run -p 389:389 -p 636:636 \
--env LDAP_ORGANISATION="BigChem" \
--env LDAP_DOMAIN="ochem.eu" \
--env LDAP_ADMIN_PASSWORD="demo" \
--volume `pwd`/test-ldap-server.ldif:/container/service/slapd/assets/config/bootstrap/ldif/custom/test-ldap-server.ldif \
--name ochem-ldap-test --rm \
osixia/openldap:1.5.0 --copy-service --loglevel debug
