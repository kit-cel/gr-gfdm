%template(preamble_generator_sptr) boost::shared_ptr<gr::gfdm::preamble_generator>;
%pythoncode %{
preamble_generator_sptr.__repr__ = lambda self: "<preamble generator>"
preamble_generator = preamble_generator.make;
%}
