# Running the Test LDAP Server

The test server uses the `osixia/openldap` Docker image that can initialize and run an easy-to-use LDAP server. You can use the `run_ldap.sh` script to run this image in a Docker container. The container is initialized with some default data about fictitious users of a fictitious organization. You can find the login credentials of the test users in the `test-ldap-server.ldif` file in this directory. When the server is running, you can test the connection by running the `test_ldap_db.sh` script and it should print out information about the organization and its users.

Once the server is running, you can configure OChem to use it for authentication of the users within the organization by setting the following variables in the configuration file:

```
ochem.login.{ldap.ochem-test.url} = ldap://localhost:389
ochem.login.{ldap.ochem-test.security_principal} = cn=admin,dc=ochem,dc=eu
ochem.login.{ldap.ochem-test.security_credentials} = demo
ochem.login.{ldap.ochem-test.root_domain} = ochem.eu
```

The `ochem-test` LDAP server should then appear in the web interface of OChem if LDAP is chosen as the login identity provider. You can then choose any user from the `test-ldap-server.ldif` file and OChem should create a user account for them the first time they successfully authenticate with the server. They can always access this account when they authenticate with the given LDAP server in the future.
